# Reproducibility notes

This repository is assembled from the actual files used for the manuscript analysis. No missing analytical results have been reconstructed from assumptions.

## Source data

The original LINCS L1000 Level 5 matrix is the GEO GSE92742 file:

`GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`

Because this is a large public source file, it is not stored in this repository and must be obtained from GEO accession GSE92742.

The following GSE92742 metadata files used by the analysis are included in `data/`:

- `GSE92742_Broad_LINCS_gene_info.txt.gz`
- `GSE92742_Broad_LINCS_pert_info.txt.gz`
- `GSE92742_Broad_LINCS_sig_info.txt.gz`

These repository copies were verified against the corresponding files used in the original analysis.

## Data paths

The repository versions of `09_extraer_level5.py` and `47_extraer_expresion_E2_temporal.py` use relative paths to the Level 5 GCTX file and the required repository data files.

For reproducible execution, the downloaded Level 5 GCTX file should be placed at:

`data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`

The GCTX file is excluded from Git because of its size.

## Hallmark gene sets

The Hallmark analyses use:

`h.all.v2026.1.Hs.symbols.gmt`

This file corresponds to the MSigDB human Hallmark gene-set collection release 2026.1.Hs and is included in `data/`.

## Provenance

The files in `data/curation/` are original curation and signature-selection files from the working analysis. Their SHA256 hashes were verified against the corresponding original analysis files during repository preparation, confirming that the repository copies are unchanged.

The GSE92742 metadata files and the Hallmark GMT file included in `data/` were likewise checked against the files used in the original analysis before inclusion in the repository.
