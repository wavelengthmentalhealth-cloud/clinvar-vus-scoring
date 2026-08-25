#!/usr/bin/env python3
"""
vus_pipeline.py

Automated Variant-of-Uncertain-Significance (VUS) triage pipeline.

Given a gene symbol, this script:
  1. Queries ClinVar (NCBI E-utilities) for every variant classified
     "Uncertain Significance" in that gene.
  2. Looks up each variant's population frequency (gnomAD), in silico
     damage predictions (CADD, SIFT, PolyPhen-2), and annotation info
     via the free myvariant.info aggregator API.
  3. Scores and ranks each variant from "leans pathogenic" to
     "leans benign" using a simple, transparent rule set.
  4. Writes a ranked CSV table and prints a summary to the console.

USAGE
-----
    pip install -r requirements.txt
    python vus_pipeline.py BRCA1
    python vus_pipeline.py MYBPC3 --max 200 --out mybpc3_vus.csv

This is a research/triage tool, not a clinical classifier. It mimics
the *logic* clinical labs use (rarity + predictor agreement + evidence)
but does not replace formal ACMG/AMP curation by a qualified lab.

NOTE: This script needs internet access to run. NCBI's E-utilities
are rate-limited to 3 requests/second without an API key (10/sec with
one) -- the script self-throttles to stay under that automatically.
"""

import argparse
import sys
import time
import csv
import xml.etree.ElementTree as ET

import requests

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MYVARIANT_API = "https://myvariant.info/v1"

# Be polite to NCBI: stay under 3 req/sec (no API key assumed).
NCBI_DELAY_SECONDS = 0.34


def ncbi_get(endpoint, params):
    """Rate-limited GET against NCBI E-utilities."""
    time.sleep(NCBI_DELAY_SECONDS)
    resp = requests.get(f"{NCBI_EUTILS}/{endpoint}", params=params, timeout=30)
    resp.raise_for_status()
    return resp


def fetch_clinvar_vus_ids(gene, max_results=100):
    """
    Step 1: Search ClinVar for all variants in `gene` classified as
    Uncertain Significance. Returns a list of ClinVar UIDs.
    """
    query = f'{gene}[gene] AND clinsig_vus[Properties] AND single_gene[Type]'
    params = {
        "db": "clinvar",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
    }
    resp = ncbi_get("esearch.fcgi", params)
    data = resp.json()
    ids = data.get("esearchresult", {}).get("idlist", [])
    print(f"[ClinVar] Found {len(ids)} VUS record(s) for {gene} "
          f"(requested up to {max_results}).")
    return ids


def fetch_clinvar_details(uids):
    """
    Step 1b: Pull structured summary details (protein change, HGVS,
    genomic location, review status) for a batch of ClinVar UIDs
    using esummary. Returns a list of dicts.
    """
    if not uids:
        return []

    variants = []
    # esummary accepts many IDs per call; keep batches modest to be safe.
    batch_size = 50
    for i in range(0, len(uids), batch_size):
        batch = uids[i:i + batch_size]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "retmode": "json",
        }
        resp = ncbi_get("esummary.fcgi", params)
        data = resp.json()
        result = data.get("result", {})
        for uid in batch:
            record = result.get(uid)
            if not record:
                continue
            variants.append(parse_clinvar_record(record))
    return variants


def parse_clinvar_record(record):
    """
    Extract the fields we care about from one ClinVar esummary record.

    NOTE: verified against a live esummary response on 2026-08-21. Two
    fields that are commonly assumed to be simple are not:
      - variation_loc[].ref / .alt are frequently blank strings, even
        though the fields exist. The reliable source for ref/alt is
        variation_set[0].canonical_spdi, which encodes them at the end
        (format: "<refseq>:<0-based pos>:<ref>:<alt>").
      - review_status is nested under germline_classification, not at
        the top level of the record.
      - protein_change is not a top-level field either; it's parsed out
        of the human-readable "title" field (e.g. "...(p.Gly72fs)"),
        since ClinVar's esummary doesn't expose it separately.
    """
    title = record.get("title", "")
    gene_info = record.get("genes", [{}])
    gene_symbol = gene_info[0].get("symbol") if gene_info else ""

    variation_set = record.get("variation_set", [{}])
    v0 = variation_set[0] if variation_set else {}
    cdna_change = v0.get("cdna_change", "")

    # Protein change isn't a separate field -- pull it out of the title,
    # e.g. "NM_000256.3(MYBPC3):c.215del (p.Gly72fs)" -> "p.Gly72fs"
    protein_change = ""
    if "(p." in title:
        protein_change = "p." + title.split("(p.")[-1].split(")")[0]

    review_status = record.get("germline_classification", {}).get("review_status", "")

    # Genomic location + ref/alt, needed to query myvariant.info.
    var_loc = v0.get("variation_loc", [])
    grch38_loc = next((loc for loc in var_loc if loc.get("assembly_name") == "GRCh38"), None)
    spdi = v0.get("canonical_spdi", "")

    hgvs_genomic = None
    variant_class = "unknown"
    if grch38_loc and spdi:
        chrom = grch38_loc.get("chr")
        start = grch38_loc.get("start")  # 1-based, matches HGVS convention
        try:
            _, _, spdi_ref, spdi_alt = spdi.split(":")
        except ValueError:
            spdi_ref, spdi_alt = "", ""

        if chrom and start and len(spdi_ref) == 1 and len(spdi_alt) == 1:
            # Simple single-base substitution -- the case our scoring
            # pipeline is built for.
            hgvs_genomic = f"chr{chrom}:g.{start}{spdi_ref}>{spdi_alt}"
            variant_class = "substitution"
        else:
            # Indel / deletion / insertion / complex variant. Building
            # correct HGVS for these needs different notation per type;
            # left unannotated here rather than guessing wrong.
            variant_class = "indel_or_complex"

    return {
        "clinvar_uid": record.get("uid", ""),
        "title": title,
        "gene": gene_symbol,
        "cdna_change": cdna_change,
        "protein_change": protein_change,
        "review_status": review_status,
        "hgvs_genomic": hgvs_genomic,
        "variant_class": variant_class,
    }


def annotate_with_myvariant(hgvs_genomic):
    """
    Step 2: Given an HGVS genomic string, query myvariant.info for
    CADD, SIFT, PolyPhen-2, and gnomAD exome allele frequency.
    Returns a dict of annotation fields (missing values as None).
    """
    empty = {
        "cadd_phred": None,
        "sift_pred": None,
        "polyphen2_pred": None,
        "gnomad_exome_af": None,
        "gnomad_genome_af": None,
    }
    if not hgvs_genomic:
        return empty  # indel/complex variant, or missing genomic location -- see variant_class field

    fields = ",".join([
        "dbnsfp.cadd.phred",
        "dbnsfp.sift.pred",
        "dbnsfp.polyphen2.hdiv.pred",
        "gnomad_exome.af.af",
        "gnomad_genome.af.af",
    ])
    try:
        resp = requests.get(
            f"{MYVARIANT_API}/variant/{hgvs_genomic}",
            # assembly=hg38 is required -- myvariant.info defaults to
            # hg19/GRCh37, but ClinVar's coordinates here are GRCh38.
            # Without this the lookup silently resolves the wrong
            # genomic position and returns no annotation.
            params={"fields": fields, "assembly": "hg38"},
            timeout=20,
        )
        if resp.status_code != 200:
            return empty
        data = resp.json()
        if isinstance(data, list):
            data = data[0] if data else {}

        dbnsfp = data.get("dbnsfp", {})
        cadd = dbnsfp.get("cadd", {})
        sift = dbnsfp.get("sift", {})
        polyphen2 = dbnsfp.get("polyphen2", {}).get("hdiv", {})

        gnomad_exome = data.get("gnomad_exome", {}).get("af", {})
        gnomad_genome = data.get("gnomad_genome", {}).get("af", {})

        def first(val):
            # dbNSFP fields are sometimes lists (multiple transcripts); take first.
            if isinstance(val, list):
                return val[0] if val else None
            return val

        return {
            "cadd_phred": first(cadd.get("phred")),
            "sift_pred": first(sift.get("pred")),
            "polyphen2_pred": first(polyphen2.get("pred")),
            "gnomad_exome_af": gnomad_exome.get("af") if isinstance(gnomad_exome, dict) else gnomad_exome,
            "gnomad_genome_af": gnomad_genome.get("af") if isinstance(gnomad_genome, dict) else gnomad_genome,
        }
    except (requests.RequestException, ValueError, IndexError):
        return empty


def score_variant(annotation):
    """
    Step 3: A simple, transparent scoring rule -- NOT a substitute for
    ACMG/AMP classification. Higher score = leans more pathogenic.

    Logic:
      +2  CADD phred >= 20 (top ~1% most damaging predicted variants)
      +1  CADD phred >= 15
      +1  SIFT predicts damaging/deleterious
      +1  PolyPhen-2 predicts damaging (probably/possibly)
      +2  Absent or ultra-rare in gnomAD (< 0.0001, i.e. <1 in 10,000)
      -2  Common in gnomAD (> 0.001, i.e. >1 in 1,000) -- too common for
          a rare dominant-disease allele in most cases
    """
    score = 0
    reasons = []

    cadd = annotation.get("cadd_phred")
    if cadd is not None:
        try:
            cadd = float(cadd)
            if cadd >= 20:
                score += 2
                reasons.append(f"CADD {cadd:.1f} (high)")
            elif cadd >= 15:
                score += 1
                reasons.append(f"CADD {cadd:.1f} (moderate)")
        except (TypeError, ValueError):
            pass

    sift = (annotation.get("sift_pred") or "").upper()
    if "D" in sift:  # SIFT uses 'D' for Damaging, 'T' for Tolerated
        score += 1
        reasons.append("SIFT: damaging")

    polyphen = (annotation.get("polyphen2_pred") or "").upper()
    if "D" in polyphen or "P" in polyphen:  # D=probably damaging, P=possibly damaging
        score += 1
        reasons.append("PolyPhen-2: damaging")

    af_values = [
        v for v in [annotation.get("gnomad_exome_af"), annotation.get("gnomad_genome_af")]
        if v is not None
    ]
    max_af = max(af_values) if af_values else None

    if max_af is None:
        score += 2
        reasons.append("Absent from gnomAD")
    elif max_af < 0.0001:
        score += 2
        reasons.append(f"Ultra-rare in gnomAD (AF={max_af:.2e})")
    elif max_af > 0.001:
        score -= 2
        reasons.append(f"Common in gnomAD (AF={max_af:.2e}) -- likely too common to be pathogenic")

    return score, "; ".join(reasons) if reasons else "Insufficient data"


def lean_label(score):
    if score >= 4:
        return "Leans pathogenic"
    elif score >= 1:
        return "Weakly leans pathogenic"
    elif score <= -2:
        return "Leans benign"
    else:
        return "Genuinely uncertain"


def generate_summary(gene, rows):
    """
    Step 4: Produce publication-ready summary statistics for a gene's
    annotated VUS set -- the numbers a Results section actually needs.
    """
    total = len(rows)
    annotated = [r for r in rows if r.get("cadd_phred") is not None
                 or r.get("gnomad_exome_af") is not None
                 or r.get("gnomad_genome_af") is not None]

    counts = {"Leans pathogenic": 0, "Weakly leans pathogenic": 0,
              "Genuinely uncertain": 0, "Leans benign": 0}
    for r in rows:
        counts[r["lean"]] = counts.get(r["lean"], 0) + 1

    cadd_values = []
    for r in rows:
        c = r.get("cadd_phred")
        if c is not None:
            try:
                cadd_values.append(float(c))
            except (TypeError, ValueError):
                pass

    absent_from_gnomad = sum(
        1 for r in rows
        if r.get("gnomad_exome_af") is None and r.get("gnomad_genome_af") is None
    )

    summary = {
        "gene": gene,
        "total_vus_pulled": total,
        "total_annotated": len(annotated),
        "pct_annotated": round(100 * len(annotated) / total, 1) if total else 0.0,
        "count_leans_pathogenic": counts.get("Leans pathogenic", 0),
        "count_weakly_pathogenic": counts.get("Weakly leans pathogenic", 0),
        "count_uncertain": counts.get("Genuinely uncertain", 0),
        "count_leans_benign": counts.get("Leans benign", 0),
        "mean_cadd": round(sum(cadd_values) / len(cadd_values), 2) if cadd_values else None,
        "n_absent_from_gnomad": absent_from_gnomad,
        "pct_absent_from_gnomad": round(100 * absent_from_gnomad / total, 1) if total else 0.0,
    }
    return summary


def write_manuscript_ready_summary(summaries, out_path="results_summary.txt"):
    """
    Writes a plain-text block of fill-in-place sentences for a manuscript
    Results section, built from real generate_summary() output. Copy the
    sentences you want directly into the paper -- every number in here
    came from an actual pipeline run, not a placeholder.
    """
    lines = []
    lines.append("=" * 70)
    lines.append("MANUSCRIPT-READY SUMMARY SENTENCES")
    lines.append("(Copy/adapt into your Results section -- every number below")
    lines.append(" is from this actual run, not a placeholder.)")
    lines.append("=" * 70)
    lines.append("")

    for s in summaries:
        lines.append(f"--- {s['gene']} ---")
        lines.append(
            f"Of {s['total_vus_pulled']} variants of uncertain significance "
            f"retrieved from ClinVar for {s['gene']}, {s['total_annotated']} "
            f"({s['pct_annotated']}%) had at least one annotation available "
            f"from gnomAD and/or dbNSFP (CADD/SIFT/PolyPhen-2) via myvariant.info."
        )
        lines.append(
            f"Applying the scoring rule described in Methods, "
            f"{s['count_leans_pathogenic']} variant(s) scored as 'leans "
            f"pathogenic,' {s['count_weakly_pathogenic']} as 'weakly leans "
            f"pathogenic,' {s['count_uncertain']} as 'genuinely uncertain,' "
            f"and {s['count_leans_benign']} as 'leans benign.'"
        )
        if s["mean_cadd"] is not None:
            lines.append(f"Mean CADD Phred score across annotated variants: {s['mean_cadd']}.")
        lines.append(
            f"{s['n_absent_from_gnomad']} variant(s) ({s['pct_absent_from_gnomad']}%) "
            f"were absent from gnomAD entirely."
        )
        lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))
    print(f"\nManuscript-ready summary written to: {out_path}")
    print("\n".join(lines))


def run_pipeline(gene, max_results=100, out_path=None):
    print(f"\n=== VUS Pipeline: {gene} ===\n")

    uids = fetch_clinvar_vus_ids(gene, max_results=max_results)
    if not uids:
        print("No VUS found. Check the gene symbol and try again.")
        return None, None

    variants = fetch_clinvar_details(uids)

    print(f"[myvariant.info] Annotating {len(variants)} variant(s) "
          f"with CADD / SIFT / PolyPhen-2 / gnomAD...")

    rows = []
    for v in variants:
        annotation = annotate_with_myvariant(v.get("hgvs_genomic"))
        score, reasons = score_variant(annotation)
        row = {**v, **annotation, "score": score, "lean": lean_label(score), "reasons": reasons}
        rows.append(row)

    rows.sort(key=lambda r: r["score"], reverse=True)

    out_path = out_path or f"{gene.upper()}_vus_ranked.csv"
    fieldnames = [
        "gene", "protein_change", "cdna_change", "clinvar_uid", "review_status",
        "cadd_phred", "sift_pred", "polyphen2_pred",
        "gnomad_exome_af", "gnomad_genome_af",
        "score", "lean", "reasons",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved ranked results to: {out_path}\n")
    print(f"{'Protein change':<20} {'Score':<6} {'Lean':<26} Reasons")
    print("-" * 100)
    for r in rows[:20]:
        pc = (r.get("protein_change") or r.get("cdna_change") or "?")[:19]
        print(f"{pc:<20} {r['score']:<6} {r['lean']:<26} {r['reasons']}")
    if len(rows) > 20:
        print(f"... and {len(rows) - 20} more (see {out_path})")

    summary = generate_summary(gene, rows)
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description="Automated ClinVar VUS triage pipeline.")
    parser.add_argument("gene", help="Gene symbol, e.g. BRCA1 or MYBPC3")
    parser.add_argument("--max", type=int, default=100, help="Max VUS to fetch (default 100)")
    parser.add_argument("--out", type=str, default=None, help="Output CSV path")
    parser.add_argument("--compare", type=str, default=None,
                         help="Second gene symbol to run in the same pass, for a "
                              "manuscript-style comparison (e.g. --compare BRCA1)")
    parser.add_argument("--summary-out", type=str, default="results_summary.txt",
                         help="Path for the manuscript-ready summary text file")
    args = parser.parse_args()

    try:
        summaries = []
        _, summary1 = run_pipeline(args.gene, max_results=args.max, out_path=args.out)
        if summary1:
            summaries.append(summary1)

        if args.compare:
            compare_out = f"{args.compare.upper()}_vus_ranked.csv"
            _, summary2 = run_pipeline(args.compare, max_results=args.max, out_path=compare_out)
            if summary2:
                summaries.append(summary2)

        if summaries:
            write_manuscript_ready_summary(summaries, out_path=args.summary_out)
    except requests.RequestException as e:
        print(f"Network error talking to a public API: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
