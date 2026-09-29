# Analysis pipeline

This document describes the organisation and execution of the analysis pipeline used for the manuscript:

**Convergent and context-dependent transcriptional programmes of vitamin D and estrogen signalling: a comparative LINCS L1000 analysis**

The analysis is implemented as 79 numbered Python scripts in the `scripts/` directory. Together, these scripts cover data exploration and curation, extraction of LINCS L1000 Level 5 expression data, construction and validation of the matched primary dataset, statistical and gene-level analyses, Hallmark gene-set analyses, context-dependent and temporal analyses, and generation of manuscript figures, tables and supplementary material.

## Working directory

The scripts are designed to be run from the **root directory of the repository**, rather than from inside the `scripts/` directory.

Several scripts generate intermediate files that are subsequently read by later scripts. These intermediate files are created in the working directory and therefore preserve the sequential dependencies of the original analysis.

## Pipeline organisation

The 79 scripts are organised into the following analytical stages:

1. **Scripts 01–08 — Data exploration, signature selection and pre-analysis curation**
2. **Scripts 09–15 — Level 5 extraction, quality control and construction of the matched primary dataset**
3. **Scripts 16–24 — Global structure, similarity and gene-level convergence analyses**
4. **Scripts 25–38 — Hallmark enrichment and Full L1000 versus Landmark robustness analyses**
5. **Scripts 39–44 — Context-dependent and complementary analyses**
6. **Scripts 45–52 — Temporal estradiol analyses**
7. **Scripts 53–58 — Main manuscript figures and tables**
8. **Scripts 59–68 — Supplementary Tables S1–S9 and their audit**
9. **Scripts 69–79 — Publication formatting and compilation of the supplementary material**

Some numbered scripts perform inspection, validation or reporting steps rather than generating a new analytical result. They are retained because they document the actual workflow used during the analysis.

## Software environment

The analysis was performed using Python 3.14.7.

The Python package versions used for the analysis are recorded in `requirements.txt`. To install these dependencies in a Python environment, run the following command from the repository root:

```bash
python -m pip install -r requirements.txt
```

For reproducibility, use the package versions specified in `requirements.txt`.

## Required input data

The analysis uses LINCS L1000 data from GEO accession **GSE92742**.

The following metadata files required by the pipeline are already included in the `data/` directory:

- `GSE92742_Broad_LINCS_gene_info.txt.gz`
- `GSE92742_Broad_LINCS_pert_info.txt.gz`
- `GSE92742_Broad_LINCS_sig_info.txt.gz`

The large Level 5 expression matrix is not included in the repository. Before running the expression-extraction steps, obtain the following file from GEO accession GSE92742:

`GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`

Place the downloaded and decompressed file at:

```text
data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
```

The `data/curation/` directory contains the curated input and signature-selection files preserved from the original analysis.

The Hallmark analyses use the MSigDB human Hallmark gene-set collection file:

`data/h.all.v2026.1.Hs.symbols.gmt`

Additional information on data provenance and file verification is provided in `docs/REPRODUCIBILITY_NOTES.md`.

## Running the pipeline

Run the scripts from the **repository root directory**. For example:

```bash
python scripts/09_extraer_level5.py
```

Do not change the working directory to `scripts/` before execution. Many stages write intermediate files to the repository root, and subsequent scripts read those files as inputs.

The numerical prefixes reflect the order in which the analytical workflow was developed and executed. The pipeline contains both analytical scripts and inspection, validation and reporting scripts; therefore, the numbered sequence also preserves quality-control and provenance steps used during the original analysis.

## Analytical workflow

### Stage 1 — Data exploration and pre-analysis curation (scripts 01–08)

These scripts explore the GSE92742 metadata, identify the perturbagens and experimental conditions of interest, restrict the candidate signatures according to the study design, and preserve the pre-expression-analysis curation and signature selection.

### Stage 2 — Level 5 extraction and primary dataset construction (scripts 09–15)

These scripts extract the selected signatures from the LINCS L1000 Level 5 matrix, perform quality-control checks, construct the matched experimental design and generate the primary 12-signature expression dataset used in subsequent analyses.

### Stage 3 — Global and gene-level analyses (scripts 16–24)

These scripts examine the global structure of the primary dataset, multivariate separation and dispersion, transcriptional similarity, gene-level programmes, directional convergence and divergence, and gene-level overlap.

### Stage 4 — Hallmark enrichment and robustness analyses (scripts 25–38)

These scripts construct gene rankings, perform Hallmark gene-set analyses, compare shared pathway-level responses and assess robustness between the complete L1000 gene space and the 978 Landmark-gene space.

### Stage 5 — Context-dependent analyses (scripts 39–44)

These scripts examine the behaviour of convergent genes and pathways across individual experimental contexts and perform complementary analyses of context heterogeneity and paired transcriptional similarity.

### Stage 6 — Temporal estradiol analysis (scripts 45–52)

These scripts construct and analyse the temporal estradiol dataset, including expression extraction, gene-order validation, temporal transcriptional analysis and Hallmark analyses.

### Stage 7 — Main manuscript outputs (scripts 53–58)

These scripts generate the five main manuscript figures and the principal manuscript tables from the analytical outputs produced in the preceding stages.

### Stage 8 — Supplementary tables and audit (scripts 59–68)

These scripts generate Supplementary Tables S1–S9 and perform consistency and integrity checks on the supplementary outputs.

### Stage 9 — Publication preparation (scripts 69–79)

These scripts prepare publication-formatted supplementary files, generate Word versions where required and compile the final supplementary material.