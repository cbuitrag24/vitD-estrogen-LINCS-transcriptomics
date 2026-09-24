# Vitamin D–estrogen LINCS transcriptomics

Reproducible LINCS L1000 analysis of convergent and context-dependent transcriptional programmes induced by estradiol, calcitriol and calcipotriol.

## Overview

This repository contains the analysis code and processed data associated with the study:

**Convergent and context-dependent transcriptional programmes of vitamin D and estrogen signalling: a comparative LINCS L1000 analysis**

The study compares transcriptional responses to canonical 17β-estradiol, calcitriol and calcipotriol using publicly available LINCS L1000 perturbational data. The analysis examines global transcriptional structure, gene-level directional consistency and overlap, Hallmark gene-set enrichment, context dependence, and temporal sensitivity.

## Data source

The analyses use publicly available LINCS L1000 Level 5 MODZ data from NCBI Gene Expression Omnibus accession **GSE92742**.

The original large Level 5 GCTX file is not redistributed in this repository. It can be obtained from the original GEO record.

## Reproducibility

This repository contains the complete analysis pipeline, comprising 79 Python scripts covering data curation, quality control, matched-design construction, statistical analyses, gene-level and pathway-level analyses, context and temporal analyses, figure generation, and preparation of manuscript and supplementary tables.

The primary analysis uses matched 24-hour perturbational signatures for estradiol, calcitriol and calcipotriol across four matched experimental contexts. Analyses are performed in the complete 12,328-gene L1000 space and, as a robustness analysis, in the 978 directly measured Landmark-gene space.

## Repository structure

The repository is organised as follows:

- `scripts/` — complete Python analysis pipeline (79 scripts), from data exploration and curation through statistical analyses, pathway analyses, figure generation and supplementary material preparation.
- `data/` — processed data files required for reproducibility that can be redistributed through the repository.
- `docs/` — supporting documentation for the analysis and reproducibility workflow.
- `requirements.txt` — Python package requirements for the analysis.
- `.gitignore` — excludes large source datasets, local environments and other files that should not be version controlled.

The original LINCS L1000 Level 5 GCTX data are not included because of their size and should be obtained directly from GEO accession GSE92742.

## Citation

Citation information will be added after publication.

## License

License information will be added before final public release of the reproducibility package.
