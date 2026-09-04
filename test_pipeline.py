"""
test_pipeline.py

Automated tests for vus_pipeline.py's pure-logic functions -- scoring,
labeling, and record parsing. These run entirely offline (no network
calls), which is what JOSS's automated-testing requirement is checking
for: that a reviewer (or CI system) can verify the tool's core logic
without needing live API access.

Run with:
    pip install pytest
    pytest tests/test_pipeline.py -v
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vus_pipeline import score_variant, lean_label, parse_clinvar_record


def test_score_high_cadd_and_absent_from_gnomad_scores_high():
    annotation = {
        "cadd_phred": 25.9,
        "sift_pred": "D",
        "polyphen2_pred": "D",
        "gnomad_exome_af": None,
        "gnomad_genome_af": None,
    }
    score, reasons = score_variant(annotation)
    assert score == 6
    assert "CADD" in reasons
    assert "Absent from gnomAD" in reasons


def test_score_common_variant_is_penalized():
    annotation = {
        "cadd_phred": 5.0,
        "sift_pred": "T",
        "polyphen2_pred": "B",
        "gnomad_exome_af": 0.05,
        "gnomad_genome_af": None,
    }
    score, reasons = score_variant(annotation)
    assert score == -2
    assert "Common in gnomAD" in reasons


def test_score_with_no_annotation_data_is_zero():
    annotation = {
        "cadd_phred": None,
        "sift_pred": None,
        "polyphen2_pred": None,
        "gnomad_exome_af": 0.0005,
        "gnomad_genome_af": None,
    }
    score, reasons = score_variant(annotation)
    assert score == 0


def test_score_moderate_cadd_only_adds_one():
    annotation = {
        "cadd_phred": 17.0,
        "sift_pred": None,
        "polyphen2_pred": None,
        "gnomad_exome_af": None,
        "gnomad_genome_af": None,
    }
    score, reasons = score_variant(annotation)
    assert score == 3


def test_score_takes_max_af_across_exome_and_genome():
    annotation = {
        "cadd_phred": None,
        "sift_pred": None,
        "polyphen2_pred": None,
        "gnomad_exome_af": 0.0000001,
        "gnomad_genome_af": 0.05,
    }
    score, reasons = score_variant(annotation)
    assert score == -2
    assert "Common in gnomAD" in reasons


def test_lean_label_boundaries():
    assert lean_label(6) == "Leans pathogenic"
    assert lean_label(4) == "Leans pathogenic"
    assert lean_label(3) == "Weakly leans pathogenic"
    assert lean_label(1) == "Weakly leans pathogenic"
    assert lean_label(0) == "Genuinely uncertain"
    assert lean_label(-1) == "Genuinely uncertain"
    assert lean_label(-2) == "Leans benign"
    assert lean_label(-5) == "Leans benign"


def _fake_clinvar_record(spdi="NC_000011.10:47342877:C:G"):
    return {
        "uid": "4874124",
        "title": "NM_000256.3(MYBPC3):c.1409G>C (p.Arg470Pro)",
        "genes": [{"symbol": "MYBPC3"}],
        "variation_set": [
            {
                "cdna_change": "c.1409G>C",
                "canonical_spdi": spdi,
                "variation_loc": [
                    {
                        "assembly_name": "GRCh38",
                        "chr": "11",
                        "start": "47342878",
                        "ref": "",
                        "alt": "",
                    },
                    {
                        "assembly_name": "GRCh37",
                        "chr": "11",
                        "start": "47364429",
                        "ref": "",
                        "alt": "",
                    },
                ],
            }
        ],
        "germline_classification": {
            "review_status": "criteria provided, single submitter",
        },
    }


def test_parse_extracts_protein_change_from_title():
    parsed = parse_clinvar_record(_fake_clinvar_record())
    assert parsed["protein_change"] == "p.Arg470Pro"


def test_parse_extracts_gene_and_review_status():
    parsed = parse_clinvar_record(_fake_clinvar_record())
    assert parsed["gene"] == "MYBPC3"
    assert parsed["review_status"] == "criteria provided, single submitter"


def test_parse_builds_hgvs_from_spdi_for_simple_substitution():
    parsed = parse_clinvar_record(_fake_clinvar_record())
    assert parsed["hgvs_genomic"] == "chr11:g.47342878C>G"
    assert parsed["variant_class"] == "substitution"


def test_parse_flags_indels_as_unsupported_rather_than_guessing():
    record = _fake_clinvar_record(spdi="NC_000011.10:47351315:CCC:CC")
    parsed = parse_clinvar_record(record)
    assert parsed["hgvs_genomic"] is None
    assert parsed["variant_class"] == "indel_or_complex"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
