import pandas as pd

# ============================================================
# Supplementary Table S6
# Leading-edge analysis of the robust convergent Hallmark core
# ============================================================

SUMMARY_FILE = "LEADING_EDGE_nucleo_robusto_resumen.csv"
SHARED_FILE = "LEADING_EDGE_genes_comunes_3.csv"
RECURRENT_FILE = "LEADING_EDGE_genes_recurrentes_nucleo.csv"

OUTPUT_A = "Supplementary_Table_S6A.csv"
OUTPUT_B = "Supplementary_Table_S6B.csv"
OUTPUT_C = "Supplementary_Table_S6C.csv"

# ------------------------------------------------------------
# Read original result files
# ------------------------------------------------------------

summary = pd.read_csv(SUMMARY_FILE)
shared = pd.read_csv(SHARED_FILE)
recurrent = pd.read_csv(RECURRENT_FILE)

# ------------------------------------------------------------
# Validate source columns
# ------------------------------------------------------------

summary_required = [
    "tratamiento",
    "espacio",
    "hallmark",
    "ES",
    "direccion",
    "pico_rank",
    "genes_set_en_ranking",
    "n_leading_edge",
]

shared_required = [
    "espacio",
    "hallmark",
    "gene",
]

recurrent_required = [
    "espacio",
    "gene",
    "n_hallmarks",
]

assert not [c for c in summary_required if c not in summary.columns]
assert not [c for c in shared_required if c not in shared.columns]
assert not [c for c in recurrent_required if c not in recurrent.columns]

# ------------------------------------------------------------
# Expected treatments, spaces and robust Hallmarks
# ------------------------------------------------------------

expected_treatments = {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}

expected_spaces = {
    "full",
    "landmark",
}

expected_hallmarks = {
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
}

# ------------------------------------------------------------
# Validate S6A source
# ------------------------------------------------------------

assert len(summary) == 18
assert set(summary["tratamiento"]) == expected_treatments
assert set(summary["espacio"]) == expected_spaces
assert set(summary["hallmark"]) == expected_hallmarks

assert (
    summary.groupby(["espacio", "tratamiento"]).size() == 3
).all()

for col in [
    "ES",
    "pico_rank",
    "genes_set_en_ranking",
    "n_leading_edge",
]:
    assert summary[col].notna().all(), f"Missing values in {col}"

# ------------------------------------------------------------
# Construct S6A
# ------------------------------------------------------------

s6a = summary[summary_required].copy()

s6a = s6a.rename(
    columns={
        "tratamiento": "Treatment",
        "espacio": "Gene_space",
        "hallmark": "Hallmark_programme",
        "direccion": "ES_direction",
        "pico_rank": "Peak_rank",
        "genes_set_en_ranking": "N_Hallmark_genes_in_ranking",
        "n_leading_edge": "N_leading_edge_genes",
    }
)

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol",
}

space_labels = {
    "full": "Full",
    "landmark": "Landmark",
}

s6a["Treatment"] = s6a["Treatment"].map(treatment_labels)
s6a["Gene_space"] = s6a["Gene_space"].map(space_labels)

space_order = {"Full": 0, "Landmark": 1}
treatment_order = {
    "Estradiol": 0,
    "Calcitriol": 1,
    "Calcipotriol": 2,
}

s6a["_space_order"] = s6a["Gene_space"].map(space_order)
s6a["_treatment_order"] = s6a["Treatment"].map(treatment_order)

s6a = (
    s6a.sort_values(
        [
            "_space_order",
            "Hallmark_programme",
            "_treatment_order",
        ],
        kind="stable",
    )
    .drop(columns=["_space_order", "_treatment_order"])
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Validate shared three-treatment leading-edge genes
# ------------------------------------------------------------

assert set(shared["espacio"]) == expected_spaces
assert set(shared["hallmark"]) == expected_hallmarks

shared_counts = (
    shared.groupby(["espacio", "hallmark"])
    .size()
    .to_dict()
)

expected_shared_counts = {
    ("full", "HALLMARK_E2F_TARGETS"): 55,
    ("full", "HALLMARK_G2M_CHECKPOINT"): 43,
    ("full", "HALLMARK_MYC_TARGETS_V1"): 30,
    ("landmark", "HALLMARK_E2F_TARGETS"): 20,
    ("landmark", "HALLMARK_G2M_CHECKPOINT"): 11,
    ("landmark", "HALLMARK_MYC_TARGETS_V1"): 9,
}

assert shared_counts == expected_shared_counts

# ------------------------------------------------------------
# Construct S6B
# ------------------------------------------------------------

s6b = shared[shared_required].copy()

s6b = s6b.rename(
    columns={
        "espacio": "Gene_space",
        "hallmark": "Hallmark_programme",
        "gene": "Gene",
    }
)

s6b["Gene_space"] = s6b["Gene_space"].map(space_labels)

s6b["_space_order"] = s6b["Gene_space"].map(space_order)

s6b = (
    s6b.sort_values(
        ["_space_order", "Hallmark_programme", "Gene"],
        kind="stable",
    )
    .drop(columns="_space_order")
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Validate recurrent shared-core genes
# ------------------------------------------------------------

assert set(recurrent["espacio"]) == expected_spaces
assert recurrent["n_hallmarks"].isin([2, 3]).all()

recurrent_counts = (
    recurrent.groupby(["espacio", "n_hallmarks"])
    .size()
    .to_dict()
)

expected_recurrent_counts = {
    ("full", 3): 4,
    ("full", 2): 29,
    ("landmark", 3): 2,
    ("landmark", 2): 11,
}

assert recurrent_counts == expected_recurrent_counts

# ------------------------------------------------------------
# Construct S6C
# ------------------------------------------------------------

s6c = recurrent[recurrent_required].copy()

s6c = s6c.rename(
    columns={
        "espacio": "Gene_space",
        "gene": "Gene",
        "n_hallmarks": "N_robust_Hallmark_programmes",
    }
)

s6c["Gene_space"] = s6c["Gene_space"].map(space_labels)

s6c["_space_order"] = s6c["Gene_space"].map(space_order)

s6c = (
    s6c.sort_values(
        [
            "_space_order",
            "N_robust_Hallmark_programmes",
            "Gene",
        ],
        ascending=[True, False, True],
        kind="stable",
    )
    .drop(columns="_space_order")
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------

assert len(s6a) == 18

assert len(s6b) == (
    55 + 43 + 30 + 20 + 11 + 9
)

assert len(s6c) == (
    4 + 29 + 2 + 11
)

# ------------------------------------------------------------
# Save without recalculating or rounding original results
# ------------------------------------------------------------

s6a.to_csv(OUTPUT_A, index=False, encoding="utf-8-sig")
s6b.to_csv(OUTPUT_B, index=False, encoding="utf-8-sig")
s6c.to_csv(OUTPUT_C, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S6 generated successfully.")

print(f"\nS6A rows: {len(s6a)}")
print("Treatment-specific leading-edge summaries.")

print(f"\nS6B rows: {len(s6b)}")
print("Three-treatment shared leading-edge genes.")

print("\nShared leading-edge genes by gene space and Hallmark:")
print(
    s6b.groupby(["Gene_space", "Hallmark_programme"])
    .size()
    .to_string()
)

print(f"\nS6C rows: {len(s6c)}")
print("Genes recurrent in at least two robust Hallmark programmes.")

print("\nRecurrent genes by gene space and number of Hallmarks:")
print(
    s6c.groupby(
        ["Gene_space", "N_robust_Hallmark_programmes"]
    )
    .size()
    .to_string()
)

print(f"\nOutput files:")
print(OUTPUT_A)
print(OUTPUT_B)
print(OUTPUT_C)