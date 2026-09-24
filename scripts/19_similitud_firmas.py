import itertools
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics.pairwise import cosine_similarity

# ============================================================
# SIMILITUD ENTRE FIRMAS TRANSCRIPCIONALES
# Dataset primario LINCS L1000
# 24 h - 10 uM
#
# Comparaciones dentro de cada bloque/celula:
#   estradiol vs calcitriol
#   estradiol vs calcipotriol
#   calcitriol vs calcipotriol
#
# Analisis:
#   1. 12,328 genes
#   2. 978 genes landmark
#
# Metricas:
#   Pearson
#   Spearman
#   similitud coseno
# ============================================================

ARCHIVO_EXPRESION = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_METADATA = "dataset_primario_metadata.csv"
ARCHIVO_GENES = "dataset_primario_genes.csv"

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

PARES = [
    ("estradiol", "calcitriol"),
    ("estradiol", "calcipotriol"),
    ("calcitriol", "calcipotriol")
]

print("=" * 72)
print("SIMILITUD ENTRE FIRMAS TRANSCRIPCIONALES")
print("=" * 72)

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------

npz = np.load(ARCHIVO_EXPRESION)

if "expression" in npz.files:
    X = npz["expression"]
elif "expresion" in npz.files:
    X = npz["expresion"]
else:
    raise ValueError(
        "No se encontro la matriz de expresion en el NPZ."
    )

meta = pd.read_csv(ARCHIVO_METADATA)
genes = pd.read_csv(ARCHIVO_GENES)

print("\nMatriz completa:")
print(X.shape)

if X.shape != (12, 12328):
    raise ValueError(
        f"Se esperaba una matriz (12, 12328) y se obtuvo {X.shape}"
    )

if len(meta) != 12:
    raise ValueError(
        "La metadata debe contener exactamente 12 firmas."
    )

if len(genes) != 12328:
    raise ValueError(
        "La metadata de genes debe contener 12328 genes."
    )

# ------------------------------------------------------------
# 2. Identificar bloque
# ------------------------------------------------------------

if "block_id" in meta.columns:
    BLOQUE = "block_id"
else:
    BLOQUE = "cell_id"

print("\nVariable de bloque:")
print(BLOQUE)

# ------------------------------------------------------------
# 3. Verificar diseño
# ------------------------------------------------------------

tabla = pd.crosstab(
    meta[BLOQUE],
    meta["pert_iname"]
)

print("\nDiseño bloque x tratamiento:")
print(tabla.to_string())

if tabla.shape != (4, 3):
    raise ValueError(
        "El diseño no corresponde a 4 bloques x 3 tratamientos."
    )

if not (tabla == 1).all().all():
    raise ValueError(
        "Debe existir exactamente una firma "
        "por tratamiento dentro de cada bloque."
    )

print("\nDiseño pareado: OK")

# ------------------------------------------------------------
# 4. Funciones de similitud
# ------------------------------------------------------------

def calcular_similitud(vector1, vector2):

    pearson = pearsonr(
        vector1,
        vector2
    ).statistic

    spearman = spearmanr(
        vector1,
        vector2
    ).statistic

    cosine = cosine_similarity(
        vector1.reshape(1, -1),
        vector2.reshape(1, -1)
    )[0, 0]

    return (
        float(pearson),
        float(spearman),
        float(cosine)
    )


# ------------------------------------------------------------
# 5. Comparaciones dentro de cada bloque
# ------------------------------------------------------------

def analizar_espacio(
    Xdatos,
    nombre_analisis
):

    resultados = []

    bloques = list(
        pd.unique(meta[BLOQUE])
    )

    for bloque in bloques:

        idx_bloque = np.where(
            meta[BLOQUE]
            .astype(str)
            .to_numpy()
            == str(bloque)
        )[0]

        for tratamiento1, tratamiento2 in PARES:

            idx1 = [
                i for i in idx_bloque
                if str(
                    meta.iloc[i]["pert_iname"]
                ) == tratamiento1
            ]

            idx2 = [
                i for i in idx_bloque
                if str(
                    meta.iloc[i]["pert_iname"]
                ) == tratamiento2
            ]

            if len(idx1) != 1 or len(idx2) != 1:
                raise ValueError(
                    f"Problema en bloque {bloque}: "
                    f"{tratamiento1} vs {tratamiento2}"
                )

            i = idx1[0]
            j = idx2[0]

            pearson, spearman, cosine = (
                calcular_similitud(
                    Xdatos[i, :],
                    Xdatos[j, :]
                )
            )

            resultados.append(
                {
                    "analysis": nombre_analisis,
                    "block": bloque,
                    "cell_id":
                        meta.iloc[i]["cell_id"],
                    "treatment_1":
                        tratamiento1,
                    "treatment_2":
                        tratamiento2,
                    "comparison":
                        tratamiento1
                        + "_vs_"
                        + tratamiento2,
                    "pearson":
                        pearson,
                    "spearman":
                        spearman,
                    "cosine":
                        cosine
                }
            )

    return pd.DataFrame(resultados)


# ------------------------------------------------------------
# 6. Espacio completo
# ------------------------------------------------------------

resultados_full = analizar_espacio(
    X,
    "full_12328_genes"
)

# ------------------------------------------------------------
# 7. Genes landmark
# ------------------------------------------------------------

genes["pr_is_lm"] = pd.to_numeric(
    genes["pr_is_lm"],
    errors="coerce"
)

mask_landmark = (
    genes["pr_is_lm"].to_numpy()
    == 1
)

X_landmark = X[
    :,
    mask_landmark
]

print("\nNumero de genes landmark:")
print(X_landmark.shape[1])

if X_landmark.shape[1] != 978:
    raise ValueError(
        "Se esperaban exactamente 978 genes landmark."
    )

resultados_landmark = analizar_espacio(
    X_landmark,
    "landmark_978_genes"
)

# ------------------------------------------------------------
# 8. Unir resultados
# ------------------------------------------------------------

resultados = pd.concat(
    [
        resultados_full,
        resultados_landmark
    ],
    ignore_index=True
)

resultados.to_csv(
    "SIMILITUD_firmas_por_bloque.csv",
    index=False
)

# ------------------------------------------------------------
# 9. Resumen por comparación
# ------------------------------------------------------------

resumen = (
    resultados
    .groupby(
        [
            "analysis",
            "comparison"
        ]
    )
    .agg(
        n_blocks=("block", "count"),

        pearson_mean=("pearson", "mean"),
        pearson_median=("pearson", "median"),
        pearson_min=("pearson", "min"),
        pearson_max=("pearson", "max"),

        spearman_mean=("spearman", "mean"),
        spearman_median=("spearman", "median"),
        spearman_min=("spearman", "min"),
        spearman_max=("spearman", "max"),

        cosine_mean=("cosine", "mean"),
        cosine_median=("cosine", "median"),
        cosine_min=("cosine", "min"),
        cosine_max=("cosine", "max")
    )
    .reset_index()
)

resumen.to_csv(
    "SIMILITUD_resumen_comparaciones.csv",
    index=False
)

# ------------------------------------------------------------
# 10. Mostrar resultados completos
# ------------------------------------------------------------

print("\n" + "=" * 72)
print("SIMILITUD POR BLOQUE")
print("=" * 72)

columnas_mostrar = [
    "analysis",
    "cell_id",
    "comparison",
    "pearson",
    "spearman",
    "cosine"
]

print(
    resultados[
        columnas_mostrar
    ].to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 11. Mostrar resumen
# ------------------------------------------------------------

print("\n" + "=" * 72)
print("RESUMEN - MEDIANAS")
print("=" * 72)

resumen_simple = resumen[
    [
        "analysis",
        "comparison",
        "pearson_median",
        "spearman_median",
        "cosine_median"
    ]
]

print(
    resumen_simple.to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 12. Comprobar cuál pareja es más similar
#     dentro de cada contexto
# ------------------------------------------------------------

ganadores = []

for analisis in resultados["analysis"].unique():

    sub_analisis = resultados[
        resultados["analysis"]
        == analisis
    ]

    for celula in sub_analisis["cell_id"].unique():

        sub = sub_analisis[
            sub_analisis["cell_id"]
            == celula
        ].copy()

        for metrica in [
            "pearson",
            "spearman",
            "cosine"
        ]:

            fila = sub.loc[
                sub[metrica].idxmax()
            ]

            ganadores.append(
                {
                    "analysis": analisis,
                    "cell_id": celula,
                    "metric": metrica,
                    "most_similar_pair":
                        fila["comparison"],
                    "similarity":
                        fila[metrica]
                }
            )

ganadores = pd.DataFrame(
    ganadores
)

ganadores.to_csv(
    "SIMILITUD_pareja_mas_similar_por_celula.csv",
    index=False
)

print("\n" + "=" * 72)
print("PAREJA MAS SIMILAR EN CADA CELULA")
print("=" * 72)

print(
    ganadores.to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 13. Concordancia entre espacio completo y landmark
# ------------------------------------------------------------

full = resultados[
    resultados["analysis"]
    == "full_12328_genes"
][
    [
        "cell_id",
        "comparison",
        "pearson",
        "spearman",
        "cosine"
    ]
].copy()

landmark = resultados[
    resultados["analysis"]
    == "landmark_978_genes"
][
    [
        "cell_id",
        "comparison",
        "pearson",
        "spearman",
        "cosine"
    ]
].copy()

comparacion_robustez = full.merge(
    landmark,
    on=[
        "cell_id",
        "comparison"
    ],
    suffixes=(
        "_full",
        "_landmark"
    )
)

comparacion_robustez.to_csv(
    "SIMILITUD_robustez_full_vs_landmark.csv",
    index=False
)

# ------------------------------------------------------------
# 14. Informe metodologico
# ------------------------------------------------------------

with open(
    "SIMILITUD_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "SIMILITUD ENTRE FIRMAS TRANSCRIPCIONALES\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "Dataset primario: 12 firmas Level 5 MODZ; "
        "4 bloques experimentales completos; "
        "3 tratamientos por bloque; "
        "24 h; 10 uM.\n\n"
    )

    f.write(
        "Las comparaciones se realizaron "
        "exclusivamente dentro de cada bloque "
        "experimental para minimizar "
        "confusion por contexto.\n\n"
    )

    f.write(
        "Se compararon tres pares:\n"
    )

    f.write(
        "1. estradiol vs calcitriol\n"
    )

    f.write(
        "2. estradiol vs calcipotriol\n"
    )

    f.write(
        "3. calcitriol vs calcipotriol\n\n"
    )

    f.write(
        "La similitud transcriptomica se "
        "cuantifico mediante correlacion "
        "de Pearson, correlacion de Spearman "
        "y similitud coseno.\n\n"
    )

    f.write(
        "Los analisis se realizaron en el "
        "espacio completo de 12,328 genes "
        "y se repitieron utilizando los "
        "978 genes landmark directamente "
        "medidos por L1000 como analisis "
        "de robustez.\n"
    )

print("\n" + "=" * 72)
print("Archivos creados:")
print(" - SIMILITUD_firmas_por_bloque.csv")
print(" - SIMILITUD_resumen_comparaciones.csv")
print(" - SIMILITUD_pareja_mas_similar_por_celula.csv")
print(" - SIMILITUD_robustez_full_vs_landmark.csv")
print(" - SIMILITUD_metodologia.txt")

print(
    "\nSIMILITUD DE FIRMAS "
    "FINALIZADA CORRECTAMENTE"
)