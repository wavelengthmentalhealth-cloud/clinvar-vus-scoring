---
title: 'ClinVar VUS Scoring Pipeline: An Automated Tool for Triaging Variants of Uncertain Significance Using Population Frequency and In Silico Prediction'
tags:
  - Python
  - genomics
  - bioinformatics
  - ClinVar
  - variant classification
  - genetic variants
authors:
  - name: Shiven Kaura
    affiliation: 1
affiliations:
  - name: North Park Secondary School
    index: 1
date: 3 September 2026
bibliography: paper.bib
---

# Summary

Clinical genetic testing routinely identifies DNA sequence changes -- variants
-- whose relationship to disease is not yet known. These are classified as
"Variants of Uncertain Significance" (VUS) and are not actionable in patient
care until reclassified. ClinVar, the NIH's public archive of clinically
observed genetic variants, hosts a large and growing number of VUS, but
resolving any individual variant currently requires a researcher to manually
cross-reference several separate public resources: ClinVar itself for the
existing classification and review history, gnomAD for population allele
frequency, and one or more in silico prediction tools (CADD, SIFT,
PolyPhen-2) for a computational estimate of functional impact. `ClinVar VUS
Scoring Pipeline` automates this workflow. Given a gene symbol, it retrieves
every ClinVar VUS for that gene via the NCBI E-utilities API, annotates each
variant with CADD, SIFT, PolyPhen-2, and gnomAD population frequency data in
a single call to the myvariant.info aggregation API [@Xin2016], applies a
transparent, pre-specified scoring rule, and outputs a ranked CSV table
along with plain-language summary statistics suitable for direct use in a
manuscript Results section. The scoring rule has been validated against
known-Pathogenic and known-Benign ClinVar variants, giving users a measured
sensitivity and specificity rather than an unverified heuristic.

# Statement of need

Triaging a VUS toward a more confident classification is a well-established
process in principle -- combine rarity in reference population data with
agreement among computational predictors -- but in practice it is done
variant-by-variant, by hand, through several different websites with
different query interfaces and coordinate systems. For a single variant this
is a manageable, if tedious, task. For systematically reviewing the dozens
to hundreds of VUS that can exist for a single clinically important gene, it
becomes impractical without automation, and impractical tasks tend not to
get done: many genes carry large VUS backlogs precisely because no
individual researcher has the time to manually triage them all.

`ClinVar VUS Scoring Pipeline` fills this gap for researchers, genetic
counselors, students, and citizen scientists who want a reproducible,
free, and fully automated first-pass triage of a gene's VUS burden,
without requiring institutional database access or paid API keys. It
deliberately does not attempt to replace formal ACMG/AMP clinical variant
classification [@Richards2015], which incorporates additional evidence types
(segregation analysis, functional assay data, case-level co-occurrence)
that are not uniformly available through free public APIs. Instead, it is
designed as a triage instrument: its ranked output is intended to direct
limited follow-up effort -- a closer literature search, a candidate for
functional validation -- toward the variants a transparent, checkable rule
identifies as most notable, across an entire gene's VUS set at once rather
than one variant at a time.

The tool was developed and validated as part of a comparative analysis of
VUS in *MYBPC3* (hypertrophic cardiomyopathy) and *BRCA1* (hereditary
breast/ovarian cancer). Critically, the scoring rule was calibrated against
known-Pathogenic and known-Benign variants in both genes, revealing that it
achieves high sensitivity (100% in this calibration set) but limited
specificity (17-30%) -- a measured, honest performance characterization
that is itself a useful contribution for anyone considering this class of
frequency-plus-prediction triage rule.

# Functionality

- Retrieves all ClinVar records classified "Uncertain Significance" for a
  given gene symbol, using the NCBI E-utilities `esearch`/`esummary`
  endpoints [@Landrum2018].
- Parses each record's genomic coordinates and canonical SPDI notation to
  construct HGVS genomic notation for single-nucleotide substitutions;
  insertions, deletions, and other complex variants are retained in the
  output but explicitly flagged as unannotated rather than given an
  incorrect or guessed annotation.
- Queries the myvariant.info API [@Xin2016] for each substitution variant to
  retrieve CADD Phred score [@Rentzsch2019], SIFT prediction
  [@NgHenikoff2003], PolyPhen-2 prediction [@Adzhubei2010], and gnomAD exome
  and genome population allele frequencies [@Karczewski2020].
- Applies a transparent, additive scoring rule combining predictor
  agreement and population rarity into a single interpretable score and
  category label.
- Includes a calibration mode that runs the identical pipeline against
  known-Pathogenic and known-Benign ClinVar records for the same gene,
  reporting sensitivity and specificity so users can evaluate the rule's
  real performance rather than trusting it blindly.
- Outputs a ranked CSV per gene and plain-text summaries suitable for
  direct use in a written report.
- Self-throttles all NCBI API requests to remain within published
  unauthenticated rate limits.
- Ships as both a command-line Python script and a Jupyter/Google Colab
  notebook, the latter requiring no local installation.
- Includes an automated test suite covering the scoring rule, label
  boundaries, and record-parsing logic against realistic fixture data, run
  entirely offline.

# Acknowledgements

The author thanks the maintainers of ClinVar, gnomAD, and myvariant.info
for providing free, well-documented public APIs without which a tool like
this would not be possible for an independent researcher to build.

# References
