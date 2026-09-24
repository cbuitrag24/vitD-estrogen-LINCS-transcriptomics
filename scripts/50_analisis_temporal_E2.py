import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, wilcoxon

print("=" * 95)
print("ANÁLISIS TEMPORAL E2 — 6 H vs 24 H")
print("7 BLOQUES EXPERIMENTALES APAREADOS")
print("=" * 95)

# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_npz = Path("dataset_E2_temporal_14x12328.npz")
archivo_meta = Path("dataset_E2_temporal_metadata.csv")
archivo_genes = Path("dataset_E2_temporal_genes_ordenados.csv")

for archivo in [archivo_npz, archivo_meta, archivo_genes]:
    if not archivo.exists():
        raise FileNotFoundError(f"No se encontró: {archivo}")

# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

datos = np.load(archivo_npz)

expresion = datos["expresion"]

meta = pd.read_csv(archivo_meta)

genes = pd.read_csv(
    archivo_genes,
    dtype={"pr_gene_id": str}
)

meta["tiempo"] = pd.to_numeric(
    meta["tiempo"],
    errors="raise"
).astype(int)

print("\nDimensión expresión:", expresion.shape)
print("Metadatos:", len(meta))
print("Genes:", len(genes))

if expresion.shape != (14, 12328):
    raise ValueError("Dimensión de expresión inesperada.")

if len(meta) != 14:
    raise ValueError("Se esperaban 14 filas de metadatos.")

if len(genes) != 12328:
    raise ValueError("Se esperaban 12328 genes.")

# ------------------------------------------------------------
# 3. Confirmar que el orden de las firmas coincide
# ------------------------------------------------------------

sig_npz = [str(x) for x in datos["sig_ids"]]
sig_meta = meta["sig_id"].astype(str).tolist()

if sig_npz != sig_meta:
    raise ValueError(
        "El orden de sig_id del NPZ no coincide con los metadatos."
    )

print("Orden de firmas verificado.")

# ------------------------------------------------------------
# 4. Identificar los 7 pares
# ------------------------------------------------------------

pares = []

for (cell, bloque), grupo in meta.groupby(
    ["cell_id", "bloque"],
    sort=True
):

    if len(grupo) != 2:
        raise ValueError(
            f"{cell} / {bloque} no contiene exactamente 2 firmas."
        )

    tiempos = sorted(grupo["tiempo"].tolist())

    if tiempos != [6, 24]:
        raise ValueError(
            f"{cell} / {bloque} no contiene exactamente 6 h y 24 h."
        )

    idx6 = grupo.index[
        grupo["tiempo"] == 6
    ][0]

    idx24 = grupo.index[
        grupo["tiempo"] == 24
    ][0]

    pares.append(
        {
            "cell_id": cell,
            "bloque": bloque,
            "idx6": idx6,
            "idx24": idx24
        }
    )

print("Pares temporales verificados:", len(pares))

# ------------------------------------------------------------
# 5. Similitud global 6 h vs 24 h dentro de cada bloque
# ------------------------------------------------------------

resultados = []

for par in pares:

    x6 = expresion[par["idx6"], :]
    x24 = expresion[par["idx24"], :]

    pearson = np.corrcoef(x6, x24)[0, 1]

    spearman = spearmanr(
        x6,
        x24
    ).statistic

    cosine = np.dot(x6, x24) / (
        np.linalg.norm(x6) *
        np.linalg.norm(x24)
    )

    diferencia_abs_media = np.mean(
        np.abs(x24 - x6)
    )

    resultados.append(
        {
            "cell_id": par["cell_id"],
            "bloque": par["bloque"],
            "pearson_6h_24h": pearson,
            "spearman_6h_24h": spearman,
            "cosine_6h_24h": cosine,
            "mean_abs_difference": diferencia_abs_media
        }
    )

res = pd.DataFrame(resultados)

print("\n" + "=" * 95)
print("SIMILITUD GLOBAL ENTRE 6 H Y 24 H")
print("=" * 95)

print(
    res.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ------------------------------------------------------------
# 6. Resumen de similitud
# ------------------------------------------------------------

print("\n" + "=" * 95)
print("RESUMEN DE SIMILITUD")
print("=" * 95)

for columna in [
    "pearson_6h_24h",
    "spearman_6h_24h",
    "cosine_6h_24h"
]:

    print(
        f"{columna}: "
        f"media={res[columna].mean():.4f}, "
        f"mediana={res[columna].median():.4f}, "
        f"mín={res[columna].min():.4f}, "
        f"máx={res[columna].max():.4f}"
    )

# ------------------------------------------------------------
# 7. Magnitud global de respuesta a 6 h y 24 h
# ------------------------------------------------------------

magnitudes = []

for par in pares:

    x6 = expresion[par["idx6"], :]
    x24 = expresion[par["idx24"], :]

    magnitudes.append(
        {
            "cell_id": par["cell_id"],
            "bloque": par["bloque"],
            "mean_abs_6h": np.mean(np.abs(x6)),
            "mean_abs_24h": np.mean(np.abs(x24))
        }
    )

mag = pd.DataFrame(magnitudes)

mag["delta_24h_minus_6h"] = (
    mag["mean_abs_24h"] -
    mag["mean_abs_6h"]
)

print("\n" + "=" * 95)
print("MAGNITUD GLOBAL DE LA RESPUESTA")
print("=" * 95)

print(
    mag.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ------------------------------------------------------------
# 8. Test apareado exacto/apropiado para n pequeño
# ------------------------------------------------------------

if np.allclose(
    mag["mean_abs_6h"],
    mag["mean_abs_24h"]
):
    estadistico = np.nan
    p_wilcoxon = 1.0
else:
    prueba = wilcoxon(
        mag["mean_abs_24h"],
        mag["mean_abs_6h"],
        alternative="two-sided",
        method="auto"
    )

    estadistico = prueba.statistic
    p_wilcoxon = prueba.pvalue

print("\nComparación apareada de magnitud global:")
print("Wilcoxon statistic:", estadistico)
print("p bilateral:", p_wilcoxon)

# ------------------------------------------------------------
# 9. Sensibilidad Landmark
# ------------------------------------------------------------

lm_mask = (
    pd.to_numeric(
        genes["pr_is_lm"],
        errors="coerce"
    )
    .fillna(0)
    .astype(int)
    .to_numpy()
    == 1
)

print("\nGenes Landmark:", lm_mask.sum())

resultados_lm = []

for par in pares:

    x6 = expresion[
        par["idx6"],
        lm_mask
    ]

    x24 = expresion[
        par["idx24"],
        lm_mask
    ]

    resultados_lm.append(
        {
            "cell_id": par["cell_id"],
            "bloque": par["bloque"],
            "pearson": np.corrcoef(x6, x24)[0, 1],
            "spearman": spearmanr(x6, x24).statistic,
            "cosine": np.dot(x6, x24) /
                      (np.linalg.norm(x6) * np.linalg.norm(x24))
        }
    )

res_lm = pd.DataFrame(resultados_lm)

print("\n" + "=" * 95)
print("SIMILITUD TEMPORAL — LANDMARK")
print("=" * 95)

print(
    res_lm.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

print("\nMedias Landmark:")
print(
    res_lm[
        ["pearson", "spearman", "cosine"]
    ].mean().to_string()
)

# ------------------------------------------------------------
# 10. Guardar
# ------------------------------------------------------------

res.to_csv(
    "E2_temporal_similitud_global.csv",
    index=False
)

mag.to_csv(
    "E2_temporal_magnitud_global.csv",
    index=False
)

res_lm.to_csv(
    "E2_temporal_similitud_landmark.csv",
    index=False
)

print("\n" + "=" * 95)
print("ARCHIVOS GENERADOS")
print("=" * 95)

print("- E2_temporal_similitud_global.csv")
print("- E2_temporal_magnitud_global.csv")
print("- E2_temporal_similitud_landmark.csv")

print("\nANÁLISIS TEMPORAL GLOBAL COMPLETADO")