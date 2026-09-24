import numpy as np
import pandas as pd

# ============================================================
# CONTROL DE CALIDAD DEL DATASET LEVEL 5
# ============================================================

print("=" * 60)
print("CONTROL DE CALIDAD - LINCS LEVEL 5")
print("=" * 60)

# ------------------------------------------------------------
# 1. Cargar matriz y metadatos
# ------------------------------------------------------------

datos = np.load("expresion_level5_213x12328.npz")

expresion = datos["expression"]
sig_ids_npz = datos["sig_ids"].astype(str)
gene_ids_npz = datos["gene_ids"].astype(str)

firmas = pd.read_csv("metadata_firmas_level5.csv")
genes = pd.read_csv(
    "metadata_genes_level5.csv",
    dtype={"pr_gene_id": str}
)

print("\nForma matriz:", expresion.shape)
print("Firmas metadata:", len(firmas))
print("Genes metadata:", len(genes))

# ------------------------------------------------------------
# 2. Comprobar correspondencia de identificadores
# ------------------------------------------------------------

if not np.array_equal(
    sig_ids_npz,
    firmas["sig_id"].astype(str).to_numpy()
):
    raise ValueError(
        "Los sig_id del NPZ y metadata_firmas_level5.csv "
        "no estan en el mismo orden."
    )

if not np.array_equal(
    gene_ids_npz,
    genes["pr_gene_id"].astype(str).to_numpy()
):
    raise ValueError(
        "Los gene_id del NPZ y metadata_genes_level5.csv "
        "no estan en el mismo orden."
    )

print("\nCorrespondencia de sig_id: OK")
print("Correspondencia de gene_id: OK")

# ------------------------------------------------------------
# 3. Valores faltantes o no finitos
# ------------------------------------------------------------

n_nan = np.isnan(expresion).sum()
n_inf = np.isinf(expresion).sum()

print("\nValores NaN:", int(n_nan))
print("Valores infinitos:", int(n_inf))

print("Minimo:", float(np.nanmin(expresion)))
print("Maximo:", float(np.nanmax(expresion)))
print("Media:", float(np.nanmean(expresion)))
print("Desvio estandar:", float(np.nanstd(expresion)))

# ------------------------------------------------------------
# 4. Genes landmark
# ------------------------------------------------------------

genes["pr_is_lm"] = pd.to_numeric(
    genes["pr_is_lm"],
    errors="coerce"
)

n_landmark = int((genes["pr_is_lm"] == 1).sum())
n_no_landmark = int((genes["pr_is_lm"] != 1).sum())

print("\nGenes landmark:", n_landmark)
print("Genes no landmark:", n_no_landmark)

# ------------------------------------------------------------
# 5. Resumen por compuesto
# ------------------------------------------------------------

print("\nFirmas por compuesto:")
print(
    firmas["pert_iname"]
    .value_counts()
    .to_string()
)

# ------------------------------------------------------------
# 6. Resumen por clase
# ------------------------------------------------------------

if "hormone_class" in firmas.columns:
    print("\nFirmas por clase:")
    print(
        firmas["hormone_class"]
        .value_counts()
        .to_string()
    )

# ------------------------------------------------------------
# 7. Resumen por linea celular y compuesto
# ------------------------------------------------------------

tabla_celula = pd.crosstab(
    firmas["cell_id"],
    firmas["pert_iname"]
)

tabla_celula.to_csv(
    "QC_firmas_por_celula_compuesto.csv"
)

print("\nNumero de lineas celulares:",
      firmas["cell_id"].nunique())

print("\nFirmas por linea celular y compuesto:")
print(tabla_celula.to_string())

# ------------------------------------------------------------
# 8. Dosis
# ------------------------------------------------------------

firmas["pert_dose"] = pd.to_numeric(
    firmas["pert_dose"],
    errors="coerce"
)

tabla_dosis = (
    firmas.groupby(
        ["pert_iname", "cell_id", "pert_dose_unit"]
    )["pert_dose"]
    .agg(["count", "nunique", "min", "max"])
    .reset_index()
)

tabla_dosis.to_csv(
    "QC_resumen_dosis.csv",
    index=False
)

# ------------------------------------------------------------
# 9. Firmas por pert_id
# ------------------------------------------------------------

tabla_pertid = (
    firmas.groupby(
        ["pert_iname", "pert_id"]
    )
    .size()
    .reset_index(name="n_firmas")
)

tabla_pertid.to_csv(
    "QC_firmas_por_pertid.csv",
    index=False
)

# ------------------------------------------------------------
# 10. Estadisticas por firma
# ------------------------------------------------------------

qc_firmas = firmas.copy()

qc_firmas["mean_z"] = np.mean(
    expresion,
    axis=1
)

qc_firmas["sd_z"] = np.std(
    expresion,
    axis=1
)

qc_firmas["median_z"] = np.median(
    expresion,
    axis=1
)

qc_firmas["max_abs_z"] = np.max(
    np.abs(expresion),
    axis=1
)

qc_firmas.to_csv(
    "QC_estadisticas_por_firma.csv",
    index=False
)

# ------------------------------------------------------------
# 11. Resumen reproducible
# ------------------------------------------------------------

with open(
    "QC_resumen_general.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write("CONTROL DE CALIDAD - LINCS LEVEL 5\n")
    f.write("=" * 50 + "\n\n")

    f.write(
        f"Forma matriz: {expresion.shape}\n"
    )

    f.write(
        f"Firmas: {len(firmas)}\n"
    )

    f.write(
        f"Genes: {len(genes)}\n"
    )

    f.write(
        f"Genes landmark: {n_landmark}\n"
    )

    f.write(
        f"Valores NaN: {int(n_nan)}\n"
    )

    f.write(
        f"Valores infinitos: {int(n_inf)}\n"
    )

    f.write(
        f"Lineas celulares: {firmas['cell_id'].nunique()}\n"
    )

print("\nArchivos QC creados:")
print(" - QC_firmas_por_celula_compuesto.csv")
print(" - QC_resumen_dosis.csv")
print(" - QC_firmas_por_pertid.csv")
print(" - QC_estadisticas_por_firma.csv")
print(" - QC_resumen_general.txt")

print("\nCONTROL DE CALIDAD FINALIZADO CORRECTAMENTE")