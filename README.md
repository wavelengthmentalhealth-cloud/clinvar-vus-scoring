# ClinVar VUS Scoring Pipeline

Automated pipeline that pulls ClinVar variants of uncertain significance (VUS),
annotates them with CADD, SIFT, PolyPhen-2, and gnomAD population frequency via
the myvariant.info API, and ranks them for functional-study triage.

Built for a comparative analysis of MYBPC3 (hypertrophic cardiomyopathy) and
BRCA1 (hereditary breast/ovarian cancer) VUS, supporting the manuscript
"Systematic Triage of ClinVar Variants of Uncertain Significance in MYBPC3 and
BRCA1 Using Population Frequency and In Silico Prediction."

## Contents
- `vus_pipeline.py` — the pipeline, runnable from the command line
- `VUS_Pipeline.ipynb` — same pipeline, runnable in Google Colab (no setup required)
- `MYBPC3_vus_ranked.csv`, `BRCA1_vus_ranked.csv` — ranked output from the analysis reported in the manuscript
- `results_summary.txt` — plain-language summary of the run

## Usage
See the notebook for a one-click run in Colab, or:
```
pip install -r requirements.txt
python vus_pipeline.py <GENE> --max 100
```
