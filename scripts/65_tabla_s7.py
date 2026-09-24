import pandas as pd

# ============================================================
# Supplementary Table S7
# Context consistency of recurrent Landmark core genes
# ============================================================

CONSISTENCY_FILE = "CONTEXTO_genes_nucleo_consistencia.csv"
CONVERGENCE_FILE = "CONTEXTO_genes_nucleo_convergencia.csv"

OUTPUT_A = "Supplementary_Table_S7A.csv"
OUTPUT_B = "Supplementary_Table_S7B.csv"

# ------------------------------------------------------------
# Read original result files
# ------------------------------------------------------------

consistency = pd.read_csv(CONSISTENCY_FILE)
convergence = pd.read_csv(CONVERGENCE_FILE)

# ------------------------------------------------------------
# Validate S7A source columns
# ------------------------------------------------------------

consistency_required = [
    "tratamiento",
    "gene",
    "n_contextos",
    "n_negativos",
    "n_positivos",
    "n_cero",
    "media_modz",
    "mediana_modz",
    "min_modz",
    "max_modz",
    "clasificacion",
]

missing = [
    c for c in consistency_required
    if c not in consistency.columns
]
assert not missing, f"Missing consistency columns: {missing}"

# ------------------------------------------------------------
# Validate S7A source structure
# ------------------------------------------------------------

assert len(consistency) == 39

expected_treatments = {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}

assert set(consistency["tratamiento"]) == expected_treatments

counts = consistency.groupby("tratamiento").size().to_dict()

assert counts == {
    "estradiol": 13,
    "calcitriol": 13,
    "calcipotriol": 13,
}

assert (consistency["n_contextos"] == 4).all()

assert (
    consistency["n_negativos"]
    + consistency["n_positivos"]
    + consistency["n_cero"]
    == consistency["n_contextos"]
).all()

# Verify the classifications directly represented in the source
classification_counts = (
    consistency.groupby(["tratamiento", "clasificacion"])
    .size()
    .to_dict()
)

expected_classification_counts = {
    ("calcipotriol", "negativo_4de4"): 4,
    ("calcipotriol", "predominio_negativo_3de4"): 9,
    ("calcitriol", "negativo_4de4"): 4,
    ("calcitriol", "predominio_negativo_3de4"): 9,
    ("estradiol", "negativo_4de4"): 1,
    ("estradiol", "predominio_negativo_3de4"): 7,
    ("estradiol", "heterogeneo"): 5,
}

assert classification_counts == expected_classification_counts

# ------------------------------------------------------------
# Construct S7A
# ------------------------------------------------------------

s7a = consistency[consistency_required].copy()

s7a = s7a.rename(
    columns={
        "tratamiento": "Treatment",
        "gene": "Gene",
        "n_contextos": "N_contexts",
        "n_negativos": "N_negative",
        "n_positivos": "N_positive",
        "n_cero": "N_zero",
        "media_modz": "Mean_MODZ",
        "mediana_modz": "Median_MODZ",
        "min_modz": "Min_MODZ",
        "max_modz": "Max_MODZ",
        "clasificacion": "Context_consistency",
    }
)

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol",
}

s7a["Treatment"] = s7a["Treatment"].map(treatment_labels)

treatment_order = {
    "Estradiol": 0,
    "Calcitriol": 1,
    "Calcipotriol": 2,
}

s7a["_treatment_order"] = s7a["Treatment"].map(treatment_order)

s7a = (
    s7a.sort_values(
        ["Gene", "_treatment_order"],
        kind="stable",
    )
    .drop(columns="_treatment_order")
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Validate S7B source columns
# ------------------------------------------------------------

convergence_required = [
    "gene",
    "direccion_media",
    "estradiol_media",
    "calcitriol_media",
    "calcipotriol_media",
    "estradiol_clasificacion",
    "calcitriol_clasificacion",
    "calcipotriol_clasificacion",
]

missing = [
    c for c in convergence_required
    if c not in convergence.columns
]
assert not missing, f"Missing convergence columns: {missing}"

assert len(convergence) == 13
assert convergence["gene"].nunique() == 13

# ------------------------------------------------------------
# Construct S7B directly from source results
# ------------------------------------------------------------

s7b = convergence[convergence_required].copy()

s7b = s7b.rename(
    columns={
        "gene": "Gene",
        "direccion_media": "Mean_direction_across_treatments",
        "estradiol_media": "Estradiol_mean_MODZ",
        "calcitriol_media": "Calcitriol_mean_MODZ",
        "calcipotriol_media": "Calcipotriol_mean_MODZ",
        "estradiol_clasificacion": "Estradiol_context_consistency",
        "calcitriol_clasificacion": "Calcitriol_context_consistency",
        "calcipotriol_clasificacion": "Calcipotriol_context_consistency",
    }
)

# ------------------------------------------------------------
# Context-robust criterion
#
# A gene is classified as context-robust when each of the
# three treatments is negative in at least 3 of 4 contexts.
#
# This flag is derived only from the existing classifications
# generated by the original context-analysis script.
# ------------------------------------------------------------

negative_consistent_classes = {
    "negativo_4de4",
    "predominio_negativo_3de4",
}

s7b["Context_robust_all_3_treatments"] = (
    s7b["Estradiol_context_consistency"].isin(
        negative_consistent_classes
    )
    & s7b["Calcitriol_context_consistency"].isin(
        negative_consistent_classes
    )
    & s7b["Calcipotriol_context_consistency"].isin(
        negative_consistent_classes
    )
)

# ------------------------------------------------------------
# Validate the context-robust subset
# ------------------------------------------------------------

expected_robust_genes = {
    "BIRC5",
    "CCNA2",
    "CCNB2",
    "CDC20",
    "CDC25A",
    "EZH2",
    "KIF2C",
    "MCM3",
}

observed_robust_genes = set(
    s7b.loc[
        s7b["Context_robust_all_3_treatments"],
        "Gene",
    ]
)

assert observed_robust_genes == expected_robust_genes
assert s7b["Context_robust_all_3_treatments"].sum() == 8

# All 13 recurrent genes had a negative mean response
# under all three treatments in the original analysis.
assert (
    s7b["Mean_direction_across_treatments"]
    == "negativa_en_los_3"
).all()

# Sort for presentation only:
# context-robust genes first, then alphabetically.
s7b = (
    s7b.sort_values(
        ["Context_robust_all_3_treatments", "Gene"],
        ascending=[False, True],
        kind="stable",
    )
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Save without rounding or recalculating MODZ values
# ------------------------------------------------------------

s7a.to_csv(OUTPUT_A, index=False, encoding="utf-8-sig")
s7b.to_csv(OUTPUT_B, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S7 generated successfully.")

print(f"\nS7A rows: {len(s7a)}")
print("Gene-by-treatment context-consistency summaries.")

print("\nContext-consistency classifications by treatment:")
print(
    s7a.groupby(["Treatment", "Context_consistency"])
    .size()
    .to_string()
)

print(f"\nS7B rows: {len(s7b)}")
print("Recurrent Landmark core genes.")

print(
    "\nGenes context-robust across all three treatments "
    "(negative in >=3/4 contexts for each treatment):"
)
print(
    ", ".join(
        s7b.loc[
            s7b["Context_robust_all_3_treatments"],
            "Gene",
        ].tolist()
    )
)

print(
    "\nNumber of context-robust genes:",
    int(s7b["Context_robust_all_3_treatments"].sum()),
)

print("\nOutput files:")
print(OUTPUT_A)
print(OUTPUT_B)