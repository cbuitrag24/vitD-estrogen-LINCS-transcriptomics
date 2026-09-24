import pandas as pd

# ============================================================
# Supplementary Table S8
# Receptor-axis transcriptional analysis
# ============================================================

COVERAGE_FILE = "RECEPTORES_cobertura_genes.csv"
SUMMARY_FILE = "RECEPTORES_resumen_tratamientos.csv"
COMPARISON_FILE = "RECEPTORES_E2_vs_VitD_por_contexto.csv"

OUTPUT_A = "Supplementary_Table_S8A.csv"
OUTPUT_B = "Supplementary_Table_S8B.csv"
OUTPUT_C = "Supplementary_Table_S8C.csv"

# ------------------------------------------------------------
# Read original result files
# ------------------------------------------------------------

coverage = pd.read_csv(COVERAGE_FILE)
summary = pd.read_csv(SUMMARY_FILE)
comparison = pd.read_csv(COMPARISON_FILE)

# ============================================================
# S8A. Gene coverage and L1000 measurement type
# ============================================================

coverage_required = [
    "gene",
    "eje",
    "presente",
    "pr_is_lm",
    "tipo_medicion",
]

missing = [c for c in coverage_required if c not in coverage.columns]
assert not missing, f"Missing coverage columns: {missing}"

assert len(coverage) == 9
assert coverage["gene"].nunique() == 9
assert coverage["presente"].all()

expected_genes = {
    "VDR",
    "RXRA",
    "RXRB",
    "RXRG",
    "CYP24A1",
    "CYP27B1",
    "ESR1",
    "ESR2",
    "GPER1",
}

assert set(coverage["gene"]) == expected_genes

# Verify Landmark status directly from the source file
landmark_genes = set(
    coverage.loc[coverage["pr_is_lm"] == 1, "gene"]
)

assert landmark_genes == {"GPER1"}
assert (coverage["pr_is_lm"] == 0).sum() == 8

s8a = coverage[coverage_required].copy()

s8a = s8a.rename(
    columns={
        "gene": "Gene",
        "eje": "Receptor_axis",
        "presente": "Present_in_L1000",
        "pr_is_lm": "Landmark_status",
        "tipo_medicion": "L1000_measurement_type",
    }
)

axis_labels = {
    "vitamina_D": "Vitamin D",
    "estrogeno": "Oestrogen",
}

measurement_labels = {
    "inferido": "Inferred",
    "landmark_directo": "Landmark_direct",
}

s8a["Receptor_axis"] = s8a["Receptor_axis"].map(axis_labels)
s8a["L1000_measurement_type"] = (
    s8a["L1000_measurement_type"].map(measurement_labels)
)

s8a["Landmark_status"] = s8a["Landmark_status"].map(
    {0: "No", 1: "Yes"}
)

axis_order = {
    "Vitamin D": 0,
    "Oestrogen": 1,
}

s8a["_axis_order"] = s8a["Receptor_axis"].map(axis_order)

s8a = (
    s8a.sort_values(
        ["_axis_order", "Gene"],
        kind="stable",
    )
    .drop(columns="_axis_order")
    .reset_index(drop=True)
)

# ============================================================
# S8B. Treatment-specific receptor-axis summaries
# ============================================================

summary_required = [
    "tratamiento",
    "eje",
    "gene",
    "media_modz",
    "mediana_modz",
    "minimo",
    "maximo",
    "n_contextos",
    "n_negativos",
    "n_positivos",
    "clasificacion",
]

missing = [c for c in summary_required if c not in summary.columns]
assert not missing, f"Missing summary columns: {missing}"

assert len(summary) == 27
assert summary["gene"].nunique() == 9

expected_treatments = {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}

assert set(summary["tratamiento"]) == expected_treatments

counts = summary.groupby("tratamiento").size().to_dict()

assert counts == {
    "estradiol": 9,
    "calcitriol": 9,
    "calcipotriol": 9,
}

assert (summary["n_contextos"] == 4).all()

# With four contexts and no separate zero-count field in this
# summary, the source counts must not exceed four.
assert (
    summary["n_negativos"] + summary["n_positivos"] <= 4
).all()

s8b = summary[summary_required].copy()

s8b = s8b.rename(
    columns={
        "tratamiento": "Treatment",
        "eje": "Receptor_axis",
        "gene": "Gene",
        "media_modz": "Mean_MODZ",
        "mediana_modz": "Median_MODZ",
        "minimo": "Min_MODZ",
        "maximo": "Max_MODZ",
        "n_contextos": "N_contexts",
        "n_negativos": "N_negative",
        "n_positivos": "N_positive",
        "clasificacion": "Context_consistency",
    }
)

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol",
}

s8b["Treatment"] = s8b["Treatment"].map(treatment_labels)
s8b["Receptor_axis"] = s8b["Receptor_axis"].map(axis_labels)

# Add measurement type by exact gene match.
# This does not alter any MODZ result.
measurement_by_gene = dict(
    zip(s8a["Gene"], s8a["L1000_measurement_type"])
)

s8b["L1000_measurement_type"] = (
    s8b["Gene"].map(measurement_by_gene)
)

assert s8b["L1000_measurement_type"].notna().all()

treatment_order = {
    "Estradiol": 0,
    "Calcitriol": 1,
    "Calcipotriol": 2,
}

s8b["_axis_order"] = s8b["Receptor_axis"].map(axis_order)
s8b["_treatment_order"] = s8b["Treatment"].map(treatment_order)

s8b = (
    s8b.sort_values(
        ["_axis_order", "Gene", "_treatment_order"],
        kind="stable",
    )
    .drop(columns=["_axis_order", "_treatment_order"])
    .reset_index(drop=True)
)

# ============================================================
# S8C. Context-specific E2 versus vitamin D comparison
# ============================================================

comparison_required = [
    "gene",
    "eje",
    "cell_id",
    "estradiol_modz",
    "calcitriol_modz",
    "calcipotriol_modz",
    "promedio_vitaminaD",
    "E2_menos_VitD",
]

missing = [
    c for c in comparison_required
    if c not in comparison.columns
]
assert not missing, f"Missing comparison columns: {missing}"

assert len(comparison) == 36
assert comparison["gene"].nunique() == 9
assert comparison["cell_id"].nunique() == 4
assert set(comparison["gene"]) == expected_genes

# Each receptor-axis gene must occur once in each of four contexts
assert (
    comparison.groupby("gene")["cell_id"].nunique() == 4
).all()

expected_cells = {
    "A549",
    "HA1E",
    "PC3",
    "VCAP",
}

assert set(comparison["cell_id"]) == expected_cells

# IMPORTANT:
# promedio_vitaminaD and E2_menos_VitD are retained exactly
# as generated by the original receptor-axis analysis.
# They are not recalculated here.

s8c = comparison[comparison_required].copy()

s8c = s8c.rename(
    columns={
        "gene": "Gene",
        "eje": "Receptor_axis",
        "cell_id": "Cell_context",
        "estradiol_modz": "Estradiol_MODZ",
        "calcitriol_modz": "Calcitriol_MODZ",
        "calcipotriol_modz": "Calcipotriol_MODZ",
        "promedio_vitaminaD": "Mean_vitamin_D_MODZ",
        "E2_menos_VitD": "Estradiol_minus_vitamin_D",
    }
)

s8c["Receptor_axis"] = s8c["Receptor_axis"].map(axis_labels)

s8c["L1000_measurement_type"] = (
    s8c["Gene"].map(measurement_by_gene)
)

assert s8c["L1000_measurement_type"].notna().all()

cell_order = {
    "A549": 0,
    "HA1E": 1,
    "PC3": 2,
    "VCAP": 3,
}

s8c["_axis_order"] = s8c["Receptor_axis"].map(axis_order)
s8c["_cell_order"] = s8c["Cell_context"].map(cell_order)

s8c = (
    s8c.sort_values(
        ["_axis_order", "Gene", "_cell_order"],
        kind="stable",
    )
    .drop(columns=["_axis_order", "_cell_order"])
    .reset_index(drop=True)
)

# ============================================================
# Final validation and save
# ============================================================

assert len(s8a) == 9
assert len(s8b) == 27
assert len(s8c) == 36

assert (
    s8a["L1000_measurement_type"] == "Landmark_direct"
).sum() == 1

assert (
    s8a["L1000_measurement_type"] == "Inferred"
).sum() == 8

s8a.to_csv(OUTPUT_A, index=False, encoding="utf-8-sig")
s8b.to_csv(OUTPUT_B, index=False, encoding="utf-8-sig")
s8c.to_csv(OUTPUT_C, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S8 generated successfully.")

print(f"\nS8A rows: {len(s8a)}")
print("Receptor-axis gene coverage and L1000 measurement type.")

print("\nMeasurement type:")
print(
    s8a["L1000_measurement_type"]
    .value_counts()
    .to_string()
)

print(f"\nS8B rows: {len(s8b)}")
print("Treatment-specific receptor-axis summaries.")

print("\nRows per treatment:")
print(
    s8b.groupby("Treatment")
    .size()
    .to_string()
)

print(f"\nS8C rows: {len(s8c)}")
print("Context-specific E2 versus vitamin D comparisons.")

print("\nContexts represented:")
print(
    s8c["Cell_context"]
    .value_counts()
    .sort_index()
    .to_string()
)

print("\nOutput files:")
print(OUTPUT_A)
print(OUTPUT_B)
print(OUTPUT_C)