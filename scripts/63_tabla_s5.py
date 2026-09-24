import pandas as pd

# ============================================================
# Supplementary Table S5
# Complete Hallmark enrichment results in the Landmark gene space
# ============================================================

INPUT_FILE = "GSEA_LANDMARK_todos_tratamientos.csv"
OUTPUT_FILE = "Supplementary_Table_S5.csv"

# Read original Landmark-space GSEA results
source = pd.read_csv(INPUT_FILE)

required_columns = [
    "tratamiento",
    "hallmark",
    "n_genes_landmark",
    "ES",
    "NES",
    "p_empirico",
    "q_BH",
]

missing = [c for c in required_columns if c not in source.columns]
assert not missing, f"Missing source columns: {missing}"

# Validate original structure
assert len(source) == 105

assert set(source["tratamiento"]) == {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}

counts = source.groupby("tratamiento").size().to_dict()

assert counts == {
    "estradiol": 35,
    "calcitriol": 35,
    "calcipotriol": 35,
}

# Each treatment must contain 35 unique evaluable Hallmarks
assert (
    source.groupby("tratamiento")["hallmark"].nunique() == 35
).all()

# Landmark analysis required at least 15 genes per Hallmark
assert (source["n_genes_landmark"] >= 15).all()

# Check numerical results
for col in [
    "n_genes_landmark",
    "ES",
    "NES",
    "p_empirico",
    "q_BH",
]:
    assert source[col].notna().all(), f"Missing values in {col}"

assert source["p_empirico"].between(0, 1).all()
assert source["q_BH"].between(0, 1).all()

# Construct S5 directly from verified source
s5 = source[required_columns].copy()

s5 = s5.rename(
    columns={
        "tratamiento": "Treatment",
        "hallmark": "Hallmark_programme",
        "n_genes_landmark": "N_Landmark_genes",
        "p_empirico": "Empirical_P",
        "q_BH": "BH_adjusted_empirical_P",
    }
)

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol",
}

s5["Treatment"] = s5["Treatment"].map(treatment_labels)

# Significance flag derived from the existing adjusted empirical P
s5["BH_significant_0.05"] = (
    s5["BH_adjusted_empirical_P"] < 0.05
)

# Sort for presentation only
treatment_order = {
    "Estradiol": 0,
    "Calcitriol": 1,
    "Calcipotriol": 2,
}

s5["_treatment_order"] = s5["Treatment"].map(treatment_order)

s5 = (
    s5.sort_values(
        ["_treatment_order", "BH_adjusted_empirical_P", "Hallmark_programme"],
        kind="stable",
    )
    .drop(columns="_treatment_order")
    .reset_index(drop=True)
)

# Final validation
assert len(s5) == 105

significant_counts = (
    s5[s5["BH_significant_0.05"]]
    .groupby("Treatment")
    .size()
    .to_dict()
)

assert significant_counts == {
    "Estradiol": 9,
    "Calcitriol": 12,
    "Calcipotriol": 10,
}

# Save without rounding original numerical values
s5.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S5 generated successfully.")
print(f"Total rows: {len(s5)}")

print("\nHallmark programmes per treatment:")
print(s5.groupby("Treatment").size().to_string())

print("\nBH-significant Hallmark programmes (adjusted empirical P < 0.05):")
print(
    s5[s5["BH_significant_0.05"]]
    .groupby("Treatment")
    .size()
    .to_string()
)

print(f"\nOutput file: {OUTPUT_FILE}")