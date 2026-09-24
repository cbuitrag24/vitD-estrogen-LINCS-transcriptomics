import pandas as pd

# ============================================================
# Supplementary Table S2
# Pairwise transcriptional similarity across complementary
# similarity metrics in Full and Landmark gene spaces
# ============================================================

DESCRIPTIVE_FILE = "SIMILITUD_resumen_comparaciones.csv"
EXACT_FILE = "SIMILITUD_prueba_pareada_exacta.csv"

OUTPUT_DESCRIPTIVE = "Supplementary_Table_S2A.csv"
OUTPUT_EXACT = "Supplementary_Table_S2B.csv"

# ------------------------------------------------------------
# Part A. Descriptive pairwise similarity
# ------------------------------------------------------------

desc = pd.read_csv(DESCRIPTIVE_FILE)

# Convert the wide source table into one row per
# gene space × treatment pair × similarity metric.
rows = []

metric_map = {
    "pearson": "Pearson",
    "spearman": "Spearman",
    "cosine": "Cosine",
}

space_map = {
    "full_12328_genes": "Full",
    "landmark_978_genes": "Landmark",
}

pair_map = {
    "calcitriol_vs_calcipotriol": "Calcitriol–calcipotriol",
    "estradiol_vs_calcitriol": "Estradiol–calcitriol",
    "estradiol_vs_calcipotriol": "Estradiol–calcipotriol",
}

for _, r in desc.iterrows():
    for metric_source, metric_label in metric_map.items():
        rows.append(
            {
                "Gene_space": space_map[r["analysis"]],
                "Treatment_pair": pair_map[r["comparison"]],
                "Metric": metric_label,
                "N_blocks": int(r["n_blocks"]),
                "Mean_similarity": r[f"{metric_source}_mean"],
                "Median_similarity": r[f"{metric_source}_median"],
                "Minimum_similarity": r[f"{metric_source}_min"],
                "Maximum_similarity": r[f"{metric_source}_max"],
            }
        )

s2a = pd.DataFrame(rows)

space_order = {"Full": 0, "Landmark": 1}
pair_order = {
    "Calcitriol–calcipotriol": 0,
    "Estradiol–calcitriol": 1,
    "Estradiol–calcipotriol": 2,
}
metric_order = {"Pearson": 0, "Spearman": 1, "Cosine": 2}

s2a["_space"] = s2a["Gene_space"].map(space_order)
s2a["_pair"] = s2a["Treatment_pair"].map(pair_order)
s2a["_metric"] = s2a["Metric"].map(metric_order)

s2a = (
    s2a.sort_values(["_space", "_metric", "_pair"], kind="stable")
    .drop(columns=["_space", "_pair", "_metric"])
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Part B. Exact paired comparisons
# ------------------------------------------------------------

exact = pd.read_csv(EXACT_FILE)

# Retain the two direct comparisons used for inference.
# The contrast against the average of the two E2-containing
# pairs was complementary/descriptive and is not included here.
exact = exact[
    ~exact["contraste"].str.contains("promedio_E2", na=False)
].copy()

exact_space_map = {
    "full": "Full",
    "landmark": "Landmark",
}

contrast_map = {
    "calcitriol_vs_calcipotriol_menos_estradiol_vs_calcitriol":
        "Calcitriol–calcipotriol vs estradiol–calcitriol",
    "calcitriol_vs_calcipotriol_menos_estradiol_vs_calcipotriol":
        "Calcitriol–calcipotriol vs estradiol–calcipotriol",
}

s2b = pd.DataFrame(
    {
        "Gene_space": exact["espacio"].map(exact_space_map),
        "Metric": exact["metrica"].map(metric_map),
        "Contrast": exact["contraste"].map(contrast_map),
        "N_blocks": exact["n_bloques"].astype(int),
        "Mean_similarity_VitD_pair": exact["media_similitud_VitD"],
        "Mean_similarity_comparator": exact["media_similitud_comparador"],
        "Mean_paired_difference": exact["diferencia_media"],
        "Median_paired_difference": exact["diferencia_mediana"],
        "Blocks_VitD_pair_higher": exact["bloques_VitD_mayor"].astype(int),
        "Blocks_VitD_pair_lower": exact["bloques_VitD_menor"].astype(int),
        "N_exact_sign_permutations": exact["n_permutaciones_exactas"].astype(int),
        "Exact_two_sided_P": exact["p_exacta_bilateral"],
    }
)

s2b["_space"] = s2b["Gene_space"].map(space_order)
s2b["_metric"] = s2b["Metric"].map(metric_order)

s2b = (
    s2b.sort_values(["_space", "_metric", "Contrast"], kind="stable")
    .drop(columns=["_space", "_metric"])
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

assert len(s2a) == 18
assert len(s2b) == 12

assert set(s2a["Gene_space"]) == {"Full", "Landmark"}
assert set(s2a["Metric"]) == {"Pearson", "Spearman", "Cosine"}
assert set(s2a["N_blocks"]) == {4}

assert set(s2b["Gene_space"]) == {"Full", "Landmark"}
assert set(s2b["Metric"]) == {"Pearson", "Spearman", "Cosine"}
assert set(s2b["N_blocks"]) == {4}
assert set(s2b["N_exact_sign_permutations"]) == {16}

# No direct paired comparison should be significant
# at the nominal two-sided P < 0.05 threshold.
assert (s2b["Exact_two_sided_P"] >= 0.05).all()

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

s2a.to_csv(OUTPUT_DESCRIPTIVE, index=False, encoding="utf-8-sig")
s2b.to_csv(OUTPUT_EXACT, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S2 generated successfully.")
print(f"S2A descriptive rows: {len(s2a)}")
print(f"S2B exact-comparison rows: {len(s2b)}")
print(f"Gene spaces: {', '.join(s2a['Gene_space'].unique())}")
print(f"Metrics: {', '.join(s2a['Metric'].unique())}")
print("Matched blocks per comparison: 4")
print("Exact sign configurations per direct comparison: 16")
print(f"Minimum exact two-sided P: {s2b['Exact_two_sided_P'].min()}")
print(f"Outputs: {OUTPUT_DESCRIPTIVE}, {OUTPUT_EXACT}\n")