import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import wilcoxon

print("=" * 100)
print("HALLMARKS TEMPORALES DE E2 — 6 H vs 24 H")
print("7 BLOQUES EXPERIMENTALES APAREADOS")
print("=" * 100)

# ------------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------------

archivo_npz = Path("dataset_E2_temporal_14x12328.npz")
archivo_meta = Path("dataset_E2_temporal_metadata.csv")
archivo_genes = Path("dataset_E2_temporal_genes_ordenados.csv")

for archivo in [archivo_npz, archivo_meta, archivo_genes]:
    if not archivo.exists():
        raise FileNotFoundError(f"No se encontró: {archivo}")

# Buscar el GMT Hallmark ya utilizado previamente
posibles_gmt = [
    Path("h.all.v2026.1.Hs.symbols.gmt"),
    Path("data/h.all.v2026.1.Hs.symbols.gmt"),
]

archivo_gmt = None

for p in posibles_gmt:
    if p.exists():
        archivo_gmt = p
        break

if archivo_gmt is None:
    raise FileNotFoundError(
        "No se encontró h.all.v2026.1.Hs.symbols.gmt"
    )

print("\nGMT utilizado:")
print(archivo_gmt)

# ------------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------------

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

if expresion.shape != (14, 12328):
    raise ValueError("La matriz no tiene dimensión 14 x 12328.")

if len(meta) != 14 or len(genes) != 12328:
    raise ValueError("Dimensiones inesperadas en metadatos o genes.")

sig_npz = [str(x) for x in datos["sig_ids"]]
sig_meta = meta["sig_id"].astype(str).tolist()

if sig_npz != sig_meta:
    raise ValueError(
        "El orden de las firmas no coincide con los metadatos."
    )

print("\nDataset verificado: 14 firmas x 12328 genes.")

# ------------------------------------------------------------------
# 3. Leer Hallmarks
# ------------------------------------------------------------------

hallmarks = {}

with open(
    archivo_gmt,
    "r",
    encoding="utf-8"
) as f:

    for linea in f:

        partes = linea.rstrip("\n").split("\t")

        if len(partes) >= 3:
            nombre = partes[0]
            miembros = set(partes[2:])
            hallmarks[nombre] = miembros

objetivos = [
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
]

faltantes = [
    h for h in objetivos
    if h not in hallmarks
]

if faltantes:
    raise ValueError(
        f"Hallmarks no encontrados en GMT: {faltantes}"
    )

print("Cinco Hallmarks objetivo encontrados.")

# ------------------------------------------------------------------
# 4. Mapear símbolos de genes
# ------------------------------------------------------------------

simbolos = (
    genes["pr_gene_symbol"]
    .astype(str)
    .to_numpy()
)

# Evitar ambigüedad por símbolos duplicados
serie_simbolos = pd.Series(simbolos)

duplicados = serie_simbolos.duplicated(
    keep=False
)

print(
    "Símbolos duplicados en gene_info:",
    int(duplicados.sum())
)

# ------------------------------------------------------------------
# 5. Construir los 7 pares
# ------------------------------------------------------------------

pares = []

for (cell, bloque), grupo in meta.groupby(
    ["cell_id", "bloque"],
    sort=True
):

    if len(grupo) != 2:
        raise ValueError(
            f"{cell}/{bloque} no tiene exactamente dos firmas."
        )

    g6 = grupo[grupo["tiempo"] == 6]
    g24 = grupo[grupo["tiempo"] == 24]

    if len(g6) != 1 or len(g24) != 1:
        raise ValueError(
            f"{cell}/{bloque} no tiene exactamente un 6 h y un 24 h."
        )

    pares.append(
        {
            "cell_id": cell,
            "bloque": bloque,
            "idx6": g6.index[0],
            "idx24": g24.index[0]
        }
    )

if len(pares) != 7:
    raise ValueError(
        f"Se esperaban 7 pares y se encontraron {len(pares)}."
    )

print("Pares temporales verificados:", len(pares))

# ------------------------------------------------------------------
# 6. Calcular score descriptivo por Hallmark
# ------------------------------------------------------------------
#
# IMPORTANTE:
# Este score NO es NES de GSEA.
#
# Es la media de los valores MODZ de los genes del Hallmark
# presentes en L1000, calculada por firma.
#
# Sirve para evaluar de forma transparente la dirección y
# magnitud del programa dentro de cada par 6 h / 24 h.
# ------------------------------------------------------------------

resultados = []

for hallmark in objetivos:

    miembros = hallmarks[hallmark]

    mask = np.array([
        simbolo in miembros
        for simbolo in simbolos
    ])

    n_genes = int(mask.sum())

    if n_genes == 0:
        raise ValueError(
            f"No hay genes disponibles para {hallmark}"
        )

    for par in pares:

        score6 = float(
            np.mean(
                expresion[
                    par["idx6"],
                    mask
                ]
            )
        )

        score24 = float(
            np.mean(
                expresion[
                    par["idx24"],
                    mask
                ]
            )
        )

        resultados.append(
            {
                "hallmark": hallmark,
                "cell_id": par["cell_id"],
                "bloque": par["bloque"],
                "n_genes": n_genes,
                "score_6h": score6,
                "score_24h": score24,
                "delta_24h_minus_6h": score24 - score6
            }
        )

res = pd.DataFrame(resultados)

# ------------------------------------------------------------------
# 7. Mostrar scores por bloque
# ------------------------------------------------------------------

print("\n" + "=" * 100)
print("SCORES POR BLOQUE")
print("Media MODZ de genes pertenecientes al Hallmark")
print("=" * 100)

for hallmark in objetivos:

    temp = res[
        res["hallmark"] == hallmark
    ]

    print("\n", hallmark)
    print("-" * 100)

    print(
        temp[
            [
                "cell_id",
                "bloque",
                "n_genes",
                "score_6h",
                "score_24h",
                "delta_24h_minus_6h"
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

# ------------------------------------------------------------------
# 8. Resumen y prueba apareada
# ------------------------------------------------------------------

resumen = []

for hallmark in objetivos:

    temp = res[
        res["hallmark"] == hallmark
    ].copy()

    x6 = temp["score_6h"].to_numpy()
    x24 = temp["score_24h"].to_numpy()

    if np.allclose(x6, x24):
        stat = np.nan
        p = 1.0
    else:
        prueba = wilcoxon(
            x24,
            x6,
            alternative="two-sided",
            method="auto"
        )

        stat = float(prueba.statistic)
        p = float(prueba.pvalue)

    n_neg6 = int((x6 < 0).sum())
    n_pos6 = int((x6 > 0).sum())

    n_neg24 = int((x24 < 0).sum())
    n_pos24 = int((x24 > 0).sum())

    resumen.append(
        {
            "hallmark": hallmark,
            "n_genes": int(temp["n_genes"].iloc[0]),
            "mean_6h": float(np.mean(x6)),
            "mean_24h": float(np.mean(x24)),
            "median_6h": float(np.median(x6)),
            "median_24h": float(np.median(x24)),
            "mean_delta_24h_minus_6h":
                float(np.mean(x24 - x6)),
            "negative_6h": n_neg6,
            "positive_6h": n_pos6,
            "negative_24h": n_neg24,
            "positive_24h": n_pos24,
            "wilcoxon_statistic": stat,
            "p_two_sided": p
        }
    )

resumen = pd.DataFrame(resumen)

# ------------------------------------------------------------------
# 9. Benjamini-Hochberg para las cinco pruebas
# ------------------------------------------------------------------

pvals = resumen["p_two_sided"].to_numpy()

orden = np.argsort(pvals)
m = len(pvals)

qvals = np.empty(m, dtype=float)

prev = 1.0

for posicion in range(m - 1, -1, -1):

    idx = orden[posicion]

    rango = posicion + 1

    q = pvals[idx] * m / rango

    prev = min(prev, q)

    qvals[idx] = min(prev, 1.0)

resumen["q_BH"] = qvals

# ------------------------------------------------------------------
# 10. Resumen final
# ------------------------------------------------------------------

print("\n" + "=" * 100)
print("RESUMEN TEMPORAL DE LOS CINCO HALLMARKS")
print("=" * 100)

columnas_resumen = [
    "hallmark",
    "n_genes",
    "mean_6h",
    "mean_24h",
    "mean_delta_24h_minus_6h",
    "negative_6h",
    "negative_24h",
    "p_two_sided",
    "q_BH"
]

print(
    resumen[
        columnas_resumen
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}"
    )
)

# ------------------------------------------------------------------
# 11. Landmark-only como control de robustez
# ------------------------------------------------------------------

lm = (
    pd.to_numeric(
        genes["pr_is_lm"],
        errors="coerce"
    )
    .fillna(0)
    .astype(int)
    .to_numpy()
    == 1
)

resultados_lm = []

for hallmark in objetivos:

    miembros = hallmarks[hallmark]

    mask = np.array([
        (simbolos[i] in miembros) and lm[i]
        for i in range(len(simbolos))
    ])

    n_genes = int(mask.sum())

    for par in pares:

        if n_genes == 0:
            score6 = np.nan
            score24 = np.nan
        else:
            score6 = float(
                np.mean(
                    expresion[
                        par["idx6"],
                        mask
                    ]
                )
            )

            score24 = float(
                np.mean(
                    expresion[
                        par["idx24"],
                        mask
                    ]
                )
            )

        resultados_lm.append(
            {
                "hallmark": hallmark,
                "cell_id": par["cell_id"],
                "bloque": par["bloque"],
                "n_landmark_genes": n_genes,
                "score_6h": score6,
                "score_24h": score24,
                "delta_24h_minus_6h": score24 - score6
            }
        )

res_lm = pd.DataFrame(resultados_lm)

print("\n" + "=" * 100)
print("COBERTURA LANDMARK DE LOS CINCO HALLMARKS")
print("=" * 100)

cobertura_lm = (
    res_lm[
        ["hallmark", "n_landmark_genes"]
    ]
    .drop_duplicates()
)

print(
    cobertura_lm.to_string(index=False)
)

# ------------------------------------------------------------------
# 12. Guardar resultados
# ------------------------------------------------------------------

res.to_csv(
    "E2_temporal_5Hallmarks_por_bloque.csv",
    index=False
)

resumen.to_csv(
    "E2_temporal_5Hallmarks_resumen.csv",
    index=False
)

res_lm.to_csv(
    "E2_temporal_5Hallmarks_landmark_por_bloque.csv",
    index=False
)

print("\n" + "=" * 100)
print("ARCHIVOS GENERADOS")
print("=" * 100)

print("- E2_temporal_5Hallmarks_por_bloque.csv")
print("- E2_temporal_5Hallmarks_resumen.csv")
print("- E2_temporal_5Hallmarks_landmark_por_bloque.csv")

print("\nANÁLISIS DE HALLMARKS TEMPORALES COMPLETADO")