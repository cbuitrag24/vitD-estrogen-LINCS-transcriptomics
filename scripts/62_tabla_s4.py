import pandas as pd

# ============================================================
# Supplementary Table S4
# Complete Hallmark enrichment results in the Full gene space
# ============================================================

INPUT_FILE = "GSEA_HALLMARK_todos_tratamientos.csv"
OUTPUT_FILE = "Supplementary_Table_S4.csv"

# ------------------------------------------------------------
# Read the original Full-space GSEA results
# ------------------------------------------------------------

source = pd.read_csv(INPUT_FILE)

required_columns = [
    "tratamiento",
    "hallmark",
    "n_genes_hallmark_en_ranking",
    "ES",
    "NES",
    "p_empirico",
    "q_BH",
]

missing = [c for c in required_columns if c not in source.columns]
assert not missing, f"Missing source columns: {missing}"

# ------------------------------------------------------------
# Validate the original result structure
# ------------------------------------------------------------

assert len(source) == 150

assert set(source["tratamiento"]) == {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}

counts = source.groupby("tratamiento").size().to_dict()

assert counts == {
    "estradiol": 50,
    "calcitriol": 50,
    "calcipotriol": 50,
}

# Each treatment should contain 50 unique Hallmark programmes
assert (
    source.groupby("tratamiento")["hallmark"].nunique() == 50
).all()

# Required numerical results must be present
for col in [
    "n_genes_hallmark_en_ranking",
    "ES",
    "NES",
    "p_empirico",
    "q_BH",
]:
    assert source[col].notna().all(), f"Missing values in {col}"

# Empirical P values and BH-adjusted P values must be valid
assert source["p_empirico"].between(0, 1).all()
assert source["q_BH"].between(0, 1).all()

# ------------------------------------------------------------
# Construct S4 directly from the verified source
# No enrichment statistic is recalculated.
# ------------------------------------------------------------

s4 = source[required_columns].copy()

s4 = s4.rename(
    columns={
        "tratamiento": "Treatment",
        "hallmark": "Hallmark_programme",
        "n_genes_hallmark_en_ranking": "N_genes_in_ranking",
        "ES": "ES",
        "NES": "NES",
        "p_empirico": "Empirical_P",
        "q_BH": "BH_adjusted_empirical_P",
    }
)

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol",
}

s4["Treatment"] = s4["Treatment"].map(treatment_labels)

# Significance indicator derived only from the existing
# BH-adjusted empirical P value.
s4["BH_significant_0.05"] = (
    s4["BH_adjusted_empirical_P"] < 0.05
)

# Sort for presentation only.
treatment_order = {
    "Estradiol": 0,
    "Calcitriol": 1,
    "Calcipotriol": 2,
}

s4["_treatment_order"] = s4["Treatment"].map(treatment_order)

s4 = (
    s4.sort_values(
        ["_treatment_order", "BH_adjusted_empirical_P", "Hallmark_programme"],
        kind="stable",
    )
    .drop(columns="_treatment_order")
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------

assert len(s4) == 150

final_counts = s4.groupby("Treatment").size().to_dict()

assert final_counts == {
    "Estradiol": 50,
    "Calcitriol": 50,
    "Calcipotriol": 50,
}

# Verify the numbers of significant Hallmarks already obtained
# from the original Full-space analysis.
significant_counts = (
    s4[s4["BH_significant_0.05"]]
    .groupby("Treatment")
    .size()
    .to_dict()
)

assert significant_counts == {
    "Estradiol": 17,
    "Calcitriol": 31,
    "Calcipotriol": 22,
}

# ------------------------------------------------------------
# Save without rounding original numerical results
# ------------------------------------------------------------

s4.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S4 generated successfully.")
print(f"Total rows: {len(s4)}")
print("\nHallmark programmes per treatment:")
print(s4.groupby("Treatment").size().to_string())

print("\nBH-significant Hallmark programmes (adjusted empirical P < 0.05):")
print(
    s4[s4["BH_significant_0.05"]]
    .groupby("Treatment")
    .size()
    .to_string()
)

print(f"\nOutput file: {OUTPUT_FILE}")