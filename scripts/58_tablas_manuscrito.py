import pandas as pd
from pathlib import Path

# ============================================================
# TABLE 1
# Matched LINCS L1000 experimental design
# ============================================================

archivo = Path("dataset_primario_metadata.csv")

df = pd.read_csv(archivo)

# Orden deseado para publicación
orden_celulas = ["A549", "HA1E", "PC3", "VCAP"]
orden_tratamientos = ["estradiol", "calcitriol", "calcipotriol"]

df["cell_id"] = pd.Categorical(
    df["cell_id"],
    categories=orden_celulas,
    ordered=True
)

df["pert_iname"] = pd.Categorical(
    df["pert_iname"],
    categories=orden_tratamientos,
    ordered=True
)

df = df.sort_values(["cell_id", "pert_iname"]).copy()

# Nombres de tratamientos para el manuscrito
nombres_tratamientos = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol"
}

# Construcción de la tabla principal
tabla1 = pd.DataFrame({
    "Cell line": df["cell_id"].astype(str),
    "Treatment": df["pert_iname"].astype(str).map(nombres_tratamientos),
    "Dose (µM)": df["pert_dose"],
    "Exposure time (h)": df["pert_time"],
    "Matched block": df["block_id"],
    "Constituent profiles (n)": df["n_replicas"]
})

# Evitar mostrar 10.0 y 24.0 si son enteros
tabla1["Dose (µM)"] = tabla1["Dose (µM)"].astype(int)
tabla1["Exposure time (h)"] = tabla1["Exposure time (h)"].astype(int)
tabla1["Constituent profiles (n)"] = tabla1["Constituent profiles (n)"].astype(int)

# Guardar versión de trabajo
salida = Path("Table_1_matched_LINCS_design.csv")
tabla1.to_csv(salida, index=False, encoding="utf-8-sig")

print("=" * 90)
print("TABLE 1. MATCHED LINCS L1000 EXPERIMENTAL DESIGN")
print("=" * 90)
print()
print(tabla1.to_string(index=False))

print("\n" + "=" * 90)
print(f"Rows: {len(tabla1)}")
print(f"Cell lines: {tabla1['Cell line'].nunique()}")
print(f"Treatments: {tabla1['Treatment'].nunique()}")
print(f"Matched blocks: {tabla1['Matched block'].nunique()}")
print(f"Total constituent profiles: {tabla1['Constituent profiles (n)'].sum()}")
print("=" * 90)

print(f"\nGenerated file: {salida}")
# ============================================================
# TABLE 2
# Global transcriptional effects and similarity
# ============================================================

print("\n\n" + "=" * 100)
print("TABLE 2. GLOBAL TRANSCRIPTIONAL EFFECTS AND SIMILARITY")
print("=" * 100)

# ------------------------------------------------------------
# A. PERMANOVA
# ------------------------------------------------------------

permanova = pd.read_csv("PERMANOVA_resultados_dataset_primario.csv")

nombres_espacio = {
    "full_12328_genes": "Full (12,328 genes)",
    "landmark_978_genes": "Landmark (978 genes)"
}

tabla2a = pd.DataFrame({
    "Gene space": permanova["analysis"].map(nombres_espacio),
    "Treatment R²": permanova["R2_treatment"],
    "Treatment partial R²": permanova["partial_R2_treatment"],
    "Treatment F": permanova["F_treatment"],
    "Exact P": permanova["p_exact_treatment"],
    "Context/block R²": permanova["R2_block_context"]
})

# ------------------------------------------------------------
# B. Multivariate dispersion
# ------------------------------------------------------------

dispersion = pd.read_csv("DISPERSION_resultados_dataset_primario.csv")

tabla2b = pd.DataFrame({
    "Gene space": dispersion["analysis"].map(nombres_espacio),
    "Dispersion F": dispersion["F_dispersion"],
    "Exact P": dispersion["p_exact_dispersion"]
})

# ------------------------------------------------------------
# C. Global transcriptional similarity
# Pearson correlation across the four matched contexts
# ------------------------------------------------------------

similitud = pd.read_csv("SIMILITUD_resumen_comparaciones.csv")

nombres_comparaciones = {
    "calcitriol_vs_calcipotriol": "Calcitriol vs calcipotriol",
    "estradiol_vs_calcitriol": "Estradiol vs calcitriol",
    "estradiol_vs_calcipotriol": "Estradiol vs calcipotriol"
}

tabla2c = pd.DataFrame({
    "Gene space": similitud["analysis"].map(nombres_espacio),
    "Treatment pair": similitud["comparison"].map(nombres_comparaciones),
    "Matched contexts (n)": similitud["n_blocks"],
    "Mean Pearson r": similitud["pearson_mean"],
    "Median Pearson r": similitud["pearson_median"],
    "Minimum Pearson r": similitud["pearson_min"],
    "Maximum Pearson r": similitud["pearson_max"]
})

# Orden de las comparaciones
orden_comparaciones = [
    "Calcitriol vs calcipotriol",
    "Estradiol vs calcitriol",
    "Estradiol vs calcipotriol"
]

tabla2c["Treatment pair"] = pd.Categorical(
    tabla2c["Treatment pair"],
    categories=orden_comparaciones,
    ordered=True
)

tabla2c = tabla2c.sort_values(
    ["Gene space", "Treatment pair"]
).reset_index(drop=True)

# ------------------------------------------------------------
# Guardar las tres secciones
# ------------------------------------------------------------

tabla2a.to_csv(
    "Table_2A_PERMANOVA.csv",
    index=False,
    encoding="utf-8-sig"
)

tabla2b.to_csv(
    "Table_2B_dispersion.csv",
    index=False,
    encoding="utf-8-sig"
)

tabla2c.to_csv(
    "Table_2C_similarity.csv",
    index=False,
    encoding="utf-8-sig"
)

# ------------------------------------------------------------
# Mostrar resultados
# ------------------------------------------------------------

print("\nA. PERMANOVA")
print("-" * 100)
print(tabla2a.to_string(index=False))

print("\nB. Multivariate dispersion")
print("-" * 100)
print(tabla2b.to_string(index=False))

print("\nC. Global transcriptional similarity")
print("-" * 100)
print(tabla2c.to_string(index=False))

print("\n" + "=" * 100)
print("Generated files:")
print("Table_2A_PERMANOVA.csv")
print("Table_2B_dispersion.csv")
print("Table_2C_similarity.csv")
print("=" * 100)
# ============================================================
# FINAL PUBLICATION VERSION — TABLE 2
# ============================================================

print("\n\n" + "=" * 105)
print("FINAL TABLE 2")
print("Global transcriptional variation, multivariate dispersion")
print("and pairwise transcriptional similarity")
print("=" * 105)

# ------------------------------------------------------------
# Section A — PERMANOVA
# ------------------------------------------------------------

final_2a = tabla2a.copy()

final_2a = final_2a.rename(columns={
    "Treatment R²": "Treatment R²",
    "Treatment partial R²": "Partial R²",
    "Treatment F": "F",
    "Exact P": "Exact P",
    "Context/block R²": "Context/block R²"
})

# ------------------------------------------------------------
# Section B — Multivariate dispersion
# ------------------------------------------------------------

final_2b = tabla2b.copy()

final_2b = final_2b.rename(columns={
    "Dispersion F": "F",
    "Exact P": "Exact P"
})

# ------------------------------------------------------------
# Section C — Pearson similarity
# ------------------------------------------------------------

final_2c = tabla2c.copy()

final_2c["Pearson r range"] = (
    final_2c["Minimum Pearson r"].map(lambda x: f"{x:.3f}") +
    " to " +
    final_2c["Maximum Pearson r"].map(lambda x: f"{x:.3f}")
)

final_2c = final_2c[
    [
        "Gene space",
        "Treatment pair",
        "Matched contexts (n)",
        "Mean Pearson r",
        "Median Pearson r",
        "Pearson r range"
    ]
].copy()

# ------------------------------------------------------------
# Round values for publication display only
# Original CSV result files remain unchanged
# ------------------------------------------------------------

for col in ["Treatment R²", "Partial R²", "F",
            "Exact P", "Context/block R²"]:
    if col in final_2a.columns:
        final_2a[col] = final_2a[col].round(3)

for col in ["F", "Exact P"]:
    if col in final_2b.columns:
        final_2b[col] = final_2b[col].round(3)

for col in ["Mean Pearson r", "Median Pearson r"]:
    final_2c[col] = final_2c[col].round(3)

# ------------------------------------------------------------
# Save publication-ready working sections
# ------------------------------------------------------------

final_2a.to_csv(
    "Table_2_final_A_PERMANOVA.csv",
    index=False,
    encoding="utf-8-sig"
)

final_2b.to_csv(
    "Table_2_final_B_dispersion.csv",
    index=False,
    encoding="utf-8-sig"
)

final_2c.to_csv(
    "Table_2_final_C_similarity.csv",
    index=False,
    encoding="utf-8-sig"
)

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

print("\nA. Global transcriptional variation (PERMANOVA)")
print("-" * 105)
print(final_2a.to_string(index=False))

print("\nB. Multivariate dispersion")
print("-" * 105)
print(final_2b.to_string(index=False))

print("\nC. Pairwise transcriptional similarity")
print("-" * 105)
print(final_2c.to_string(index=False))

print("\n" + "=" * 105)
print("TABLE NOTE")
print("=" * 105)

nota = (
    "PERMANOVA was based on Euclidean distances and exact treatment "
    "permutations within matched experimental blocks (1,296 permutations). "
    "Context/block R² is reported descriptively because cellular context "
    "and experimental block were confounded in the matched design; therefore, "
    "no inferential P value was assigned to this component. Pairwise "
    "transcriptional similarity is summarised using Pearson correlation "
    "across the four matched contexts. Although calcitriol and calcipotriol "
    "showed the highest mean similarity, differences relative to "
    "estradiol-containing pairs did not reach statistical significance "
    "in exact paired tests (all two-sided P >= 0.125 in the Landmark "
    "analysis and P >= 0.250 in the full gene space). Complete Pearson, "
    "Spearman and cosine similarities and exact paired tests are provided "
    "in Supplementary Table S2."
)

print(nota)

print("\nGenerated final Table 2 files:")
print("Table_2_final_A_PERMANOVA.csv")
print("Table_2_final_B_dispersion.csv")
print("Table_2_final_C_similarity.csv")
# ============================================================
# TABLE 3A — DIRECTION-SPECIFIC GENE-LEVEL CONVERGENCE
# ============================================================

overlap = pd.read_csv(
    "ENRIQUECIMIENTO_solapamiento_resultados.csv"
)

space_names = {
    "full_12328_genes": "Full (12,328 genes)",
    "landmark_978_genes": "Landmark (978 genes)"
}

comparison_names = {
    "calcitriol_vs_calcipotriol": "Calcitriol vs calcipotriol",
    "estradiol_vs_calcitriol": "Estradiol vs calcitriol",
    "estradiol_vs_calcipotriol": "Estradiol vs calcipotriol"
}

direction_names = {
    "up": "Up",
    "down": "Down"
}

table3a = overlap[
    [
        "analysis_space",
        "comparison",
        "direction",
        "observed_overlap",
        "expected_overlap",
        "fold_enrichment",
        "q_BH"
    ]
].copy()

table3a["analysis_space"] = table3a["analysis_space"].map(space_names)
table3a["comparison"] = table3a["comparison"].map(comparison_names)
table3a["direction"] = table3a["direction"].map(direction_names)

table3a = table3a.rename(columns={
    "analysis_space": "Gene space",
    "comparison": "Treatment pair",
    "direction": "Direction",
    "observed_overlap": "Observed overlap",
    "expected_overlap": "Expected overlap",
    "fold_enrichment": "Fold enrichment",
    "q_BH": "BH-adjusted P"
})

# Publication display formatting only
table3a["Expected overlap"] = table3a["Expected overlap"].round(2)
table3a["Fold enrichment"] = table3a["Fold enrichment"].round(2)

table3a.to_csv(
    "Table_3A_gene_level_convergence.csv",
    index=False,
    encoding="utf-8-sig"
)

print("\n\n" + "=" * 110)
print("TABLE 3A. DIRECTION-SPECIFIC GENE-LEVEL CONVERGENCE")
print("=" * 110)
print(table3a.to_string(index=False))
print("=" * 110)
# ============================================================
# TABLE 3B — ROBUST HALLMARK-LEVEL CONVERGENCE
# ============================================================

hallmark_core = pd.read_csv(
    "Figure_4_Hallmark_core_values.csv"
)

treatment_names = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol"
}

hallmark_names = {
    "HALLMARK_E2F_TARGETS": "E2F targets",
    "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
    "HALLMARK_MYC_TARGETS_V1": "MYC targets V1"
}

table3b = hallmark_core[
    [
        "tratamiento",
        "hallmark",
        "NES_full",
        "q_full",
        "NES_landmark",
        "q_landmark"
    ]
].copy()

table3b["tratamiento"] = table3b["tratamiento"].map(treatment_names)
table3b["hallmark"] = table3b["hallmark"].map(hallmark_names)

table3b = table3b.rename(columns={
    "tratamiento": "Treatment",
    "hallmark": "Hallmark programme",
    "NES_full": "Full NES",
    "q_full": "Full BH-adjusted empirical P",
"NES_landmark": "Landmark NES",
"q_landmark": "Landmark BH-adjusted empirical P"
})

# Order for publication
hallmark_order = [
    "E2F targets",
    "G2M checkpoint",
    "MYC targets V1"
]

treatment_order = [
    "Estradiol",
    "Calcitriol",
    "Calcipotriol"
]

table3b["Hallmark programme"] = pd.Categorical(
    table3b["Hallmark programme"],
    categories=hallmark_order,
    ordered=True
)

table3b["Treatment"] = pd.Categorical(
    table3b["Treatment"],
    categories=treatment_order,
    ordered=True
)

table3b = table3b.sort_values(
    ["Hallmark programme", "Treatment"]
).reset_index(drop=True)

# Round NES only for publication display.
# Adjusted P values retain their original precision.
table3b["Full NES"] = table3b["Full NES"].round(3)
table3b["Landmark NES"] = table3b["Landmark NES"].round(3)

table3b.to_csv(
    "Table_3B_Hallmark_convergence.csv",
    index=False,
    encoding="utf-8-sig"
)

print("\n\n" + "=" * 110)
print("TABLE 3B. ROBUST HALLMARK-LEVEL CONVERGENCE")
print("=" * 110)
print(table3b.to_string(index=False))
print("=" * 110)