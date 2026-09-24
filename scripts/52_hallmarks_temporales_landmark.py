import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

print("=" * 100)
print("SENSIBILIDAD TEMPORAL E2 — HALLMARKS LANDMARK")
print("6 H vs 24 H — 7 BLOQUES EXPERIMENTALES APAREADOS")
print("=" * 100)

# ------------------------------------------------------------------
# 1. Cargar resultados Landmark generados por el script 51
# ------------------------------------------------------------------

archivo = "E2_temporal_5Hallmarks_landmark_por_bloque.csv"

df = pd.read_csv(archivo)

columnas_necesarias = [
    "hallmark",
    "cell_id",
    "bloque",
    "n_landmark_genes",
    "score_6h",
    "score_24h",
    "delta_24h_minus_6h"
]

faltantes = [
    c for c in columnas_necesarias
    if c not in df.columns
]

if faltantes:
    raise ValueError(
        f"Faltan columnas necesarias: {faltantes}"
    )

print("\nArchivo Landmark cargado correctamente.")
print("Filas:", len(df))

hallmarks = [
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
]

# ------------------------------------------------------------------
# 2. Verificaciones
# ------------------------------------------------------------------

if len(df) != 35:
    raise ValueError(
        f"Se esperaban 35 filas (5 Hallmarks x 7 bloques), "
        f"pero hay {len(df)}."
    )

for hallmark in hallmarks:

    temp = df[df["hallmark"] == hallmark]

    if len(temp) != 7:
        raise ValueError(
            f"{hallmark} no tiene exactamente 7 bloques."
        )

print("Verificación correcta: 5 Hallmarks x 7 bloques.")

# ------------------------------------------------------------------
# 3. Mostrar resultados Landmark por bloque
# ------------------------------------------------------------------

print("\n" + "=" * 100)
print("SCORES LANDMARK POR BLOQUE")
print("=" * 100)

for hallmark in hallmarks:

    temp = df[
        df["hallmark"] == hallmark
    ].copy()

    print("\n" + hallmark)
    print("-" * 100)

    print(
        temp[
            [
                "cell_id",
                "bloque",
                "n_landmark_genes",
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
# 4. Resumen apareado Landmark
# ------------------------------------------------------------------

resumen = []

for hallmark in hallmarks:

    temp = df[
        df["hallmark"] == hallmark
    ].copy()

    x6 = temp["score_6h"].to_numpy(dtype=float)
    x24 = temp["score_24h"].to_numpy(dtype=float)

    delta = x24 - x6

    if np.allclose(delta, 0):
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

    resumen.append(
        {
            "hallmark": hallmark,
            "n_landmark_genes":
                int(temp["n_landmark_genes"].iloc[0]),

            "mean_6h":
                float(np.mean(x6)),

            "mean_24h":
                float(np.mean(x24)),

            "median_6h":
                float(np.median(x6)),

            "median_24h":
                float(np.median(x24)),

            "mean_delta_24h_minus_6h":
                float(np.mean(delta)),

            "median_delta_24h_minus_6h":
                float(np.median(delta)),

            "delta_negative":
                int((delta < 0).sum()),

            "delta_positive":
                int((delta > 0).sum()),

            "negative_6h":
                int((x6 < 0).sum()),

            "negative_24h":
                int((x24 < 0).sum()),

            "wilcoxon_statistic":
                stat,

            "p_two_sided":
                p
        }
    )

resumen = pd.DataFrame(resumen)

# ------------------------------------------------------------------
# 5. Corrección Benjamini-Hochberg
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
# 6. Resumen final
# ------------------------------------------------------------------

print("\n" + "=" * 100)
print("RESUMEN TEMPORAL LANDMARK")
print("=" * 100)

columnas_mostrar = [
    "hallmark",
    "n_landmark_genes",
    "mean_6h",
    "mean_24h",
    "mean_delta_24h_minus_6h",
    "delta_negative",
    "delta_positive",
    "negative_6h",
    "negative_24h",
    "p_two_sided",
    "q_BH"
]

print(
    resumen[
        columnas_mostrar
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}"
    )
)

# ------------------------------------------------------------------
# 7. Resumen por célula
# ------------------------------------------------------------------
#
# Esto NO se usa como una nueva prueba estadística.
# Sólo permite comprobar que los bloques repetidos de PC3
# no oculten la dirección general por línea celular.
# ------------------------------------------------------------------

por_celula = (
    df.groupby(
        ["hallmark", "cell_id"],
        as_index=False
    )
    .agg(
        n_bloques=("bloque", "count"),
        mean_score_6h=("score_6h", "mean"),
        mean_score_24h=("score_24h", "mean"),
        mean_delta=("delta_24h_minus_6h", "mean")
    )
)

print("\n" + "=" * 100)
print("COMPROBACIÓN DESCRIPTIVA POR LÍNEA CELULAR")
print("A549 = 2 bloques; PC3 = 4 bloques; VCAP = 1 bloque")
print("=" * 100)

for hallmark in hallmarks:

    temp = por_celula[
        por_celula["hallmark"] == hallmark
    ]

    print("\n" + hallmark)
    print("-" * 100)

    print(
        temp.to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

# ------------------------------------------------------------------
# 8. Guardar
# ------------------------------------------------------------------

resumen.to_csv(
    "E2_temporal_5Hallmarks_landmark_resumen.csv",
    index=False
)

por_celula.to_csv(
    "E2_temporal_5Hallmarks_landmark_por_celula.csv",
    index=False
)

print("\n" + "=" * 100)
print("ARCHIVOS GENERADOS")
print("=" * 100)

print(
    "- E2_temporal_5Hallmarks_landmark_resumen.csv"
)

print(
    "- E2_temporal_5Hallmarks_landmark_por_celula.csv"
)

print("\nCONTROL TEMPORAL LANDMARK COMPLETADO")