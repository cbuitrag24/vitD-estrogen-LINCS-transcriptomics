import pandas as pd

# ============================================================
# Supplementary Table S9
# Temporal analysis of estradiol response: 6 h versus 24 h
# ============================================================

AUDIT_FILE = "E2_temporal_auditoria_bloques.csv"
SIMILARITY_FILE = "E2_temporal_similitud_global.csv"
MAGNITUDE_FILE = "E2_temporal_magnitud_global.csv"
FULL_HALLMARK_FILE = "E2_temporal_5Hallmarks_resumen.csv"
LANDMARK_HALLMARK_FILE = "E2_temporal_5Hallmarks_landmark_resumen.csv"

OUTPUT_A = "Supplementary_Table_S9A.csv"
OUTPUT_B = "Supplementary_Table_S9B.csv"
OUTPUT_C = "Supplementary_Table_S9C.csv"

# ============================================================
# S9A. Audit of matched 6 h / 24 h experimental blocks
# ============================================================

audit = pd.read_csv(AUDIT_FILE)

audit_required = [
    "cell_id",
    "bloque",
    "6",
    "24",
    "tiene_6h",
    "tiene_24h",
    "pareado_6h_24h",
]

missing = [c for c in audit_required if c not in audit.columns]
assert not missing, f"Missing audit columns: {missing}"

# The source audit contains 17 candidate cell/block combinations.
assert len(audit) == 17

# Exactly seven matched blocks were identified.
matched = audit[audit["pareado_6h_24h"] == True].copy()

assert len(matched) == 7

# Verify distribution across cell lines directly from source.
matched_counts = (
    matched.groupby("cell_id")
    .size()
    .to_dict()
)

assert matched_counts == {
    "A549": 2,
    "PC3": 4,
    "VCAP": 1,
}

assert set(matched["cell_id"]) == {
    "A549",
    "PC3",
    "VCAP",
}

# Construct S9A using only matched blocks.
s9a = matched[audit_required].copy()

s9a = s9a.rename(
    columns={
        "cell_id": "Cell_line",
        "bloque": "Experimental_block",
        "6": "N_6h_profiles",
        "24": "N_24h_profiles",
        "tiene_6h": "Has_6h",
        "tiene_24h": "Has_24h",
        "pareado_6h_24h": "Matched_6h_24h",
    }
)

cell_order = {
    "A549": 0,
    "PC3": 1,
    "VCAP": 2,
}

s9a["_cell_order"] = s9a["Cell_line"].map(cell_order)

s9a = (
    s9a.sort_values(
        ["_cell_order", "Experimental_block"],
        kind="stable",
    )
    .drop(columns="_cell_order")
    .reset_index(drop=True)
)

assert len(s9a) == 7

# ============================================================
# S9B. Global temporal similarity and magnitude
# ============================================================

similarity = pd.read_csv(SIMILARITY_FILE)
magnitude = pd.read_csv(MAGNITUDE_FILE)

similarity_required = [
    "cell_id",
    "bloque",
    "pearson_6h_24h",
    "spearman_6h_24h",
    "cosine_6h_24h",
    "mean_abs_difference",
]

magnitude_required = [
    "cell_id",
    "bloque",
    "mean_abs_6h",
    "mean_abs_24h",
    "delta_24h_minus_6h",
]

missing = [
    c for c in similarity_required
    if c not in similarity.columns
]
assert not missing, f"Missing similarity columns: {missing}"

missing = [
    c for c in magnitude_required
    if c not in magnitude.columns
]
assert not missing, f"Missing magnitude columns: {missing}"

# Both files must contain the same seven matched blocks.
assert len(similarity) == 7
assert len(magnitude) == 7

similarity_keys = set(
    zip(similarity["cell_id"], similarity["bloque"])
)

magnitude_keys = set(
    zip(magnitude["cell_id"], magnitude["bloque"])
)

audit_keys = set(
    zip(matched["cell_id"], matched["bloque"])
)

assert similarity_keys == audit_keys
assert magnitude_keys == audit_keys

# Merge the already calculated quantities.
s9b = similarity[similarity_required].copy()

s9b = s9b.merge(
    magnitude[magnitude_required],
    on=["cell_id", "bloque"],
    how="left",
    validate="one_to_one",
)

assert len(s9b) == 7
assert s9b["mean_abs_6h"].notna().all()
assert s9b["mean_abs_24h"].notna().all()
assert s9b["delta_24h_minus_6h"].notna().all()

s9b = s9b.rename(
    columns={
        "cell_id": "Cell_line",
        "bloque": "Experimental_block",
        "pearson_6h_24h": "Pearson_6h_vs_24h",
        "spearman_6h_24h": "Spearman_6h_vs_24h",
        "cosine_6h_24h": "Cosine_6h_vs_24h",
        "mean_abs_difference": "Mean_absolute_difference",
        "mean_abs_6h": "Mean_absolute_MODZ_6h",
        "mean_abs_24h": "Mean_absolute_MODZ_24h",
        "delta_24h_minus_6h": "Delta_24h_minus_6h",
    }
)

s9b["_cell_order"] = s9b["Cell_line"].map(cell_order)

s9b = (
    s9b.sort_values(
        ["_cell_order", "Experimental_block"],
        kind="stable",
    )
    .drop(columns="_cell_order")
    .reset_index(drop=True)
)

# ============================================================
# S9C. Temporal Hallmark analysis
# ============================================================

full = pd.read_csv(FULL_HALLMARK_FILE)
landmark = pd.read_csv(LANDMARK_HALLMARK_FILE)

full_required = [
    "hallmark",
    "n_genes",
    "mean_6h",
    "mean_24h",
    "median_6h",
    "median_24h",
    "mean_delta_24h_minus_6h",
    "negative_6h",
    "positive_6h",
    "negative_24h",
    "positive_24h",
    "wilcoxon_statistic",
    "p_two_sided",
    "q_BH",
]

landmark_required = [
    "hallmark",
    "n_landmark_genes",
    "mean_6h",
    "mean_24h",
    "median_6h",
    "median_24h",
    "mean_delta_24h_minus_6h",
    "median_delta_24h_minus_6h",
    "delta_negative",
    "delta_positive",
    "negative_6h",
    "negative_24h",
    "wilcoxon_statistic",
    "p_two_sided",
    "q_BH",
]

missing = [c for c in full_required if c not in full.columns]
assert not missing, f"Missing Full Hallmark columns: {missing}"

missing = [
    c for c in landmark_required
    if c not in landmark.columns
]
assert not missing, f"Missing Landmark Hallmark columns: {missing}"

# Exactly five Hallmark programmes in each space.
assert len(full) == 5
assert len(landmark) == 5

expected_hallmarks = {
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
}

assert set(full["hallmark"]) == expected_hallmarks
assert set(landmark["hallmark"]) == expected_hallmarks

# ------------------------------------------------------------
# Full component
# ------------------------------------------------------------

full_out = full[full_required].copy()

full_out["Gene_space"] = "Full"

full_out = full_out.rename(
    columns={
        "hallmark": "Hallmark_programme",
        "n_genes": "N_genes",
        "mean_6h": "Mean_6h",
        "mean_24h": "Mean_24h",
        "median_6h": "Median_6h",
        "median_24h": "Median_24h",
        "mean_delta_24h_minus_6h": "Mean_delta_24h_minus_6h",
        "negative_6h": "N_negative_6h",
        "positive_6h": "N_positive_6h",
        "negative_24h": "N_negative_24h",
        "positive_24h": "N_positive_24h",
        "wilcoxon_statistic": "Wilcoxon_statistic",
        "p_two_sided": "P_two_sided",
        "q_BH": "Q_BH",
    }
)

# ------------------------------------------------------------
# Landmark component
# ------------------------------------------------------------

landmark_out = landmark[landmark_required].copy()

landmark_out["Gene_space"] = "Landmark"

landmark_out = landmark_out.rename(
    columns={
        "hallmark": "Hallmark_programme",
        "n_landmark_genes": "N_Landmark_genes",
        "mean_6h": "Mean_6h",
        "mean_24h": "Mean_24h",
        "median_6h": "Median_6h",
        "median_24h": "Median_24h",
        "mean_delta_24h_minus_6h": "Mean_delta_24h_minus_6h",
        "median_delta_24h_minus_6h": "Median_delta_24h_minus_6h",
        "delta_negative": "N_delta_negative",
        "delta_positive": "N_delta_positive",
        "negative_6h": "N_negative_6h",
        "negative_24h": "N_negative_24h",
        "wilcoxon_statistic": "Wilcoxon_statistic",
        "p_two_sided": "P_two_sided",
        "q_BH": "Q_BH",
    }
)

# ------------------------------------------------------------
# Keep Full and Landmark information explicit.
# Because their source columns differ, S9C uses a common
# set of columns with NA where the source does not contain
# that quantity.
# ------------------------------------------------------------

common_columns = [
    "Gene_space",
    "Hallmark_programme",
    "N_genes",
    "N_Landmark_genes",
    "Mean_6h",
    "Mean_24h",
    "Median_6h",
    "Median_24h",
    "Mean_delta_24h_minus_6h",
    "Median_delta_24h_minus_6h",
    "N_delta_negative",
    "N_delta_positive",
    "N_negative_6h",
    "N_positive_6h",
    "N_negative_24h",
    "N_positive_24h",
    "Wilcoxon_statistic",
    "P_two_sided",
    "Q_BH",
]

for col in common_columns:
    if col not in full_out.columns:
        full_out[col] = pd.NA

    if col not in landmark_out.columns:
        landmark_out[col] = pd.NA

s9c = pd.concat(
    [
        full_out[common_columns],
        landmark_out[common_columns],
    ],
    ignore_index=True,
)

# ------------------------------------------------------------
# Validate numerical P/q ranges
# ------------------------------------------------------------

assert s9c["P_two_sided"].between(0, 1).all()
assert s9c["Q_BH"].between(0, 1).all()

# Five Hallmarks per gene space.
assert (
    s9c.groupby("Gene_space").size().to_dict()
    == {
        "Full": 5,
        "Landmark": 5,
    }
)

# ============================================================
# Save
# ============================================================

s9a.to_csv(
    OUTPUT_A,
    index=False,
    encoding="utf-8-sig",
)

s9b.to_csv(
    OUTPUT_B,
    index=False,
    encoding="utf-8-sig",
)

s9c.to_csv(
    OUTPUT_C,
    index=False,
    encoding="utf-8-sig",
)

# ============================================================
# Final report
# ============================================================

print("\nSupplementary Table S9 generated successfully.")

print(f"\nS9A rows: {len(s9a)}")
print("Matched 6 h / 24 h experimental blocks.")

print("\nMatched blocks by cell line:")
print(
    s9a.groupby("Cell_line")
    .size()
    .to_string()
)

print(f"\nS9B rows: {len(s9b)}")
print("Global temporal similarity and magnitude.")

print(f"\nS9C rows: {len(s9c)}")
print("Five Hallmark programmes in Full and Landmark spaces.")

print("\nHallmark rows by gene space:")
print(
    s9c.groupby("Gene_space")
    .size()
    .to_string()
)

print("\nOutput files:")
print(OUTPUT_A)
print(OUTPUT_B)
print(OUTPUT_C)