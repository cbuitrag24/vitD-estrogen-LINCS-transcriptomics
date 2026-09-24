import pandas as pd

# ============================================================
# Supplementary Table S3
# Directionally consistent genes across all four matched
# experimental contexts
# ============================================================

INPUT_FILE = "PROGRAMAS_genes_consenso_4o4.csv"
OUTPUT_FILE = "Supplementary_Table_S3.csv"

# Leer directamente el archivo de resultados ya generado
source = pd.read_csv(INPUT_FILE)

# ------------------------------------------------------------
# Validaciones de la fuente
# ------------------------------------------------------------

required_columns = [
    "pr_gene_id",
    "pr_gene_symbol",
    "pr_gene_title",
    "pr_is_lm",
    "treatment",
    "n_up",
    "n_down",
    "n_zero",
    "direction",
    "consistency",
    "mean_modz",
    "median_modz",
    "mean_abs_modz",
]

missing = [c for c in required_columns if c not in source.columns]
assert not missing, f"Missing source columns: {missing}"

# Este archivo debe contener exclusivamente genes 4/4.
assert set(source["consistency"].astype(str)) == {"4_of_4"}

# Tratamientos y direcciones esperados del análisis primario.
assert set(source["treatment"]) == {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}
assert set(source["direction"]) == {"up", "down"}

# Un gen 4/4 Up debe tener 4 Up, 0 Down y 0 Zero.
up = source["direction"].eq("up")
assert (source.loc[up, "n_up"] == 4).all()
assert (source.loc[up, "n_down"] == 0).all()
assert (source.loc[up, "n_zero"] == 0).all()

# Un gen 4/4 Down debe tener 4 Down, 0 Up y 0 Zero.
down = source["direction"].eq("down")
assert (source.loc[down, "n_down"] == 4).all()
assert (source.loc[down, "n_up"] == 0).all()
assert (source.loc[down, "n_zero"] == 0).all()

# Comprobar contra los conteos previamente verificados
observed_counts = (
    source.groupby(["treatment", "direction"])
    .size()
    .to_dict()
)

verified_counts = {
    ("estradiol", "up"): 782,
    ("estradiol", "down"): 402,
    ("calcitriol", "up"): 548,
    ("calcitriol", "down"): 716,
    ("calcipotriol", "up"): 1063,
    ("calcipotriol", "down"): 897,
}

assert observed_counts == verified_counts

# ------------------------------------------------------------
# Construcción de Supplementary Table S3
# No se recalcula ningún valor.
# ------------------------------------------------------------

s3 = source[required_columns].copy()

s3 = s3.rename(
    columns={
        "pr_gene_id": "Gene_ID",
        "pr_gene_symbol": "Gene_symbol",
        "pr_gene_title": "Gene_title",
        "pr_is_lm": "Landmark_status",
        "treatment": "Treatment",
        "n_up": "N_contexts_up",
        "n_down": "N_contexts_down",
        "n_zero": "N_contexts_zero",
        "direction": "Direction",
        "consistency": "Consistency",
        "mean_modz": "Mean_MODZ",
        "median_modz": "Median_MODZ",
        "mean_abs_modz": "Mean_absolute_MODZ",
    }
)

# Etiquetas para presentación
treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol",
}

direction_labels = {
    "up": "Up",
    "down": "Down",
}

s3["Treatment"] = s3["Treatment"].map(treatment_labels)
s3["Direction"] = s3["Direction"].map(direction_labels)

# Ordenar solamente para presentación.
# No altera los datos.
treatment_order = {
    "Estradiol": 0,
    "Calcitriol": 1,
    "Calcipotriol": 2,
}
direction_order = {
    "Up": 0,
    "Down": 1,
}

s3["_treatment_order"] = s3["Treatment"].map(treatment_order)
s3["_direction_order"] = s3["Direction"].map(direction_order)

s3 = (
    s3.sort_values(
        [
            "_treatment_order",
            "_direction_order",
            "Gene_symbol",
        ],
        kind="stable",
    )
    .drop(columns=["_treatment_order", "_direction_order"])
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# Validaciones finales
# ------------------------------------------------------------

assert len(s3) == len(source)
assert len(s3) == 4408

assert s3["Gene_ID"].notna().all()
assert s3["Gene_symbol"].notna().all()

# Guardar sin redondear los valores numéricos originales
s3.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S3 generated successfully.")
print(f"Total rows: {len(s3)}")
print("\nVerified 4/4 gene counts:")

summary = (
    s3.groupby(["Treatment", "Direction"])
    .size()
    .reset_index(name="N_genes")
)

print(summary.to_string(index=False))
print(f"\nOutput file: {OUTPUT_FILE}")