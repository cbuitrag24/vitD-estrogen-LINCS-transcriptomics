# Reproducibility notes

This repository is being assembled from the actual files used for the manuscript analysis. No missing analytical results are reconstructed from assumptions.

## Source data

The original LINCS L1000 Level 5 matrix is the GEO GSE92742 file:

`GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`

Because this is a large public source file, it is not stored in this repository.

The corresponding GSE92742 signature, perturbagen and gene metadata files are also public source data.

## Script audit

The complete analysis pipeline contains 79 Python scripts numbered 01–79.

Two working scripts contain a machine-specific absolute path to the Level 5 GCTX file:

- `09_extraer_level5.py`
- `47_extraer_expresion_E2_temporal.py`

For the reproducibility release, only these path declarations will be made configurable. Analytical calculations will not be changed silently.

## Hallmark gene sets

The Hallmark analyses use `h.all.v2026.1.Hs.symbols.gmt` (MSigDB human release 2026.1.Hs), as specified by the analysis scripts. Redistribution versus download instructions will be resolved before public release.

## Provenance

The files in `data/curation/` are original curation and signature-selection files from the working analysis. Their SHA256 hashes were verified against the corresponding original analysis files during repository preparation, confirming that the repository copies are unchanged.
