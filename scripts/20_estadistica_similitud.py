import itertools
import numpy as np
import pandas as pd

# ============================================================
# ESTADISTICA DE SIMILITUD TRANSCRIPCIONAL
#
# Diseño:
# 4 contextos celulares apareados
#
# Comparaciones:
# A = calcitriol vs calcipotriol
# B = estradiol vs calcitriol
# C = estradiol vs calcipotriol
#
# Spearman = métrica principal
# Pearson y coseno = análisis de robustez
#
# Debido al pequeño número de contextos (n=4),
# se utilizan permutaciones exactas apareadas.
# ============================================================

ARCHIVO = "SIMILITUD_firmas_por_bloque.csv"

P_VITD = "calcitriol_vs_calcipotriol"
P_E2_CAL = "estradiol_vs_calcitriol"
P_E2_CALI = "estradiol_vs_calcipotriol"

METRICAS = [
    "spearman",
    "pearson",
    "cosine"
]

ANALISIS = [
    "full_12328_genes",
    "landmark_978_genes"
]

print("=" * 76)
print("ESTADISTICA DE SIMILITUD TRANSCRIPCIONAL")
print("=" * 76)

datos = pd.read_csv(ARCHIVO)

# ------------------------------------------------------------
# 1. Verificaciones
# ------------------------------------------------------------

columnas_necesarias = [
    "analysis",
    "cell_id",
    "comparison",
    "spearman",
    "pearson",
    "cosine"
]

faltantes = [
    c for c in columnas_necesarias
    if c not in datos.columns
]

if faltantes:
    raise ValueError(
        f"Faltan columnas: {faltantes}"
    )

print("\nFilas cargadas:")
print(len(datos))

# ------------------------------------------------------------
# 2. Permutación exacta apareada
#
# H0: no existe diferencia sistemática entre
# las dos comparaciones.
#
# Para cada contexto, bajo H0 podemos intercambiar
# las dos observaciones.
#
# Con 4 contextos existen 2^4 = 16 asignaciones.
#
# Estadístico: diferencia media apareada.
#
# Se informa:
# - p bilateral exacto
# - p unilateral exacto para la hipótesis
#   VitD-VitD > comparación E2-VitD
# ------------------------------------------------------------

def permutacion_exacta_pareada(x, y):

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) != len(y):
        raise ValueError(
            "Los vectores apareados tienen "
            "distinta longitud."
        )

    d = x - y

    observado = np.mean(d)

    estadisticos = []

    for signos in itertools.product(
        [-1, 1],
        repeat=len(d)
    ):

        signos = np.asarray(
            signos,
            dtype=float
        )

        estadisticos.append(
            np.mean(
                signos * d
            )
        )

    estadisticos = np.asarray(
        estadisticos
    )

    # Bilateral
    p_two_sided = np.mean(
        np.abs(estadisticos)
        >=
        (abs(observado) - 1e-12)
    )

    # Unilateral:
    # x tiene mayor similitud que y
    p_greater = np.mean(
        estadisticos
        >=
        (observado - 1e-12)
    )

    return {
        "mean_difference": observado,
        "median_difference":
            np.median(d),
        "min_difference":
            np.min(d),
        "max_difference":
            np.max(d),
        "n_positive":
            int(np.sum(d > 0)),
        "n_negative":
            int(np.sum(d < 0)),
        "n_zero":
            int(np.sum(d == 0)),
        "p_exact_two_sided":
            p_two_sided,
        "p_exact_greater":
            p_greater,
        "n_permutations":
            len(estadisticos)
    }


# ------------------------------------------------------------
# 3. Construir matrices por contexto
# ------------------------------------------------------------

resultados = []
detalles = []

for analisis in ANALISIS:

    sub = datos[
        datos["analysis"]
        == analisis
    ].copy()

    print("\n" + "=" * 76)
    print(analisis)
    print("=" * 76)

    for metrica in METRICAS:

        tabla = sub.pivot(
            index="cell_id",
            columns="comparison",
            values=metrica
        )

        necesarias = [
            P_VITD,
            P_E2_CAL,
            P_E2_CALI
        ]

        for columna in necesarias:
            if columna not in tabla.columns:
                raise ValueError(
                    f"Falta {columna} "
                    f"en {analisis}, {metrica}"
                )

        tabla = tabla[
            necesarias
        ].dropna()

        if len(tabla) != 4:
            raise ValueError(
                f"Se esperaban 4 contextos "
                f"y se encontraron {len(tabla)}."
            )

        # ------------------------------------
        # Contraste 1
        # VitD-VitD vs E2-calcitriol
        # ------------------------------------

        r1 = permutacion_exacta_pareada(
            tabla[P_VITD].values,
            tabla[P_E2_CAL].values
        )

        resultados.append(
            {
                "analysis": analisis,
                "metric": metrica,
                "contrast":
                    "VitD_VitD_vs_E2_calcitriol",
                **r1
            }
        )

        # ------------------------------------
        # Contraste 2
        # VitD-VitD vs E2-calcipotriol
        # ------------------------------------

        r2 = permutacion_exacta_pareada(
            tabla[P_VITD].values,
            tabla[P_E2_CALI].values
        )

        resultados.append(
            {
                "analysis": analisis,
                "metric": metrica,
                "contrast":
                    "VitD_VitD_vs_E2_calcipotriol",
                **r2
            }
        )

        # ------------------------------------
        # Guardar diferencias contexto a contexto
        # ------------------------------------

        for celula in tabla.index:

            detalles.append(
                {
                    "analysis": analisis,
                    "metric": metrica,
                    "cell_id": celula,

                    "vitd_vitd":
                        tabla.loc[
                            celula,
                            P_VITD
                        ],

                    "e2_calcitriol":
                        tabla.loc[
                            celula,
                            P_E2_CAL
                        ],

                    "e2_calcipotriol":
                        tabla.loc[
                            celula,
                            P_E2_CALI
                        ],

                    "delta_vs_e2_calcitriol":
                        tabla.loc[
                            celula,
                            P_VITD
                        ]
                        -
                        tabla.loc[
                            celula,
                            P_E2_CAL
                        ],

                    "delta_vs_e2_calcipotriol":
                        tabla.loc[
                            celula,
                            P_VITD
                        ]
                        -
                        tabla.loc[
                            celula,
                            P_E2_CALI
                        ]
                }
            )


# ------------------------------------------------------------
# 4. Resultados
# ------------------------------------------------------------

resultados = pd.DataFrame(
    resultados
)

detalles = pd.DataFrame(
    detalles
)

resultados.to_csv(
    "ESTADISTICA_similitud_resultados.csv",
    index=False
)

detalles.to_csv(
    "ESTADISTICA_similitud_por_contexto.csv",
    index=False
)

# ------------------------------------------------------------
# 5. Mostrar Spearman
# ------------------------------------------------------------

principal = resultados[
    resultados["metric"]
    == "spearman"
].copy()

print("\n")
print("=" * 76)
print("RESULTADO PRINCIPAL: SPEARMAN")
print("=" * 76)

columnas = [
    "analysis",
    "contrast",
    "mean_difference",
    "median_difference",
    "n_positive",
    "n_negative",
    "p_exact_two_sided",
    "p_exact_greater"
]

print(
    principal[
        columnas
    ].to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 6. Mostrar diferencias por contexto
#    para landmark/Spearman
# ------------------------------------------------------------

contexto_principal = detalles[
    (
        detalles["analysis"]
        == "landmark_978_genes"
    )
    &
    (
        detalles["metric"]
        == "spearman"
    )
]

print("\n")
print("=" * 76)
print("LANDMARK + SPEARMAN: DIFERENCIAS POR CONTEXTO")
print("=" * 76)

print(
    contexto_principal[
        [
            "cell_id",
            "vitd_vitd",
            "e2_calcitriol",
            "e2_calcipotriol",
            "delta_vs_e2_calcitriol",
            "delta_vs_e2_calcipotriol"
        ]
    ].round(6).to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 7. Robustez
# ------------------------------------------------------------

print("\n")
print("=" * 76)
print("ROBUSTEZ: TODAS LAS METRICAS")
print("=" * 76)

print(
    resultados[
        [
            "analysis",
            "metric",
            "contrast",
            "mean_difference",
            "n_positive",
            "p_exact_two_sided",
            "p_exact_greater"
        ]
    ].round(6).to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 8. Nota sobre resolución del test
# ------------------------------------------------------------

print("\n")
print("=" * 76)
print("NOTA ESTADISTICA")
print("=" * 76)

print(
    "Con n=4 contextos existen solamente "
    "2^4 = 16 permutaciones exactas."
)

print(
    "Por lo tanto, la resolucion minima "
    "posible del p exacto es 1/16 = 0.0625."
)

print(
    "Los tamaños y la consistencia de los "
    "efectos entre contextos son especialmente "
    "importantes para la interpretacion."
)

# ------------------------------------------------------------
# 9. Informe metodológico
# ------------------------------------------------------------

with open(
    "ESTADISTICA_similitud_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "ESTADISTICA DE SIMILITUD "
        "TRANSCRIPCIONAL\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "Las similitudes entre pares de "
        "tratamientos se compararon dentro "
        "de los cuatro contextos celulares "
        "del diseño primario pareado.\n\n"
    )

    f.write(
        "Spearman se definió como métrica "
        "principal. Pearson y similitud "
        "coseno se utilizaron como análisis "
        "de robustez.\n\n"
    )

    f.write(
        "Los contrastes principales fueron:\n"
    )

    f.write(
        "1. calcitriol-calcipotriol frente "
        "a estradiol-calcitriol.\n"
    )

    f.write(
        "2. calcitriol-calcipotriol frente "
        "a estradiol-calcipotriol.\n\n"
    )

    f.write(
        "La inferencia utilizó una prueba "
        "de permutación exacta apareada "
        "basada en la diferencia media "
        "de similitud entre comparaciones.\n"
    )

    f.write(
        "Con cuatro contextos independientes "
        "existen 2^4 = 16 asignaciones "
        "exactas de signos.\n\n"
    )

    f.write(
        "Debido a esta limitada resolución "
        "inferencial, se informaron además "
        "las diferencias de similitud y "
        "su consistencia entre contextos.\n"
    )

print("\nArchivos creados:")
print(" - ESTADISTICA_similitud_resultados.csv")
print(" - ESTADISTICA_similitud_por_contexto.csv")
print(" - ESTADISTICA_similitud_metodologia.txt")

print(
    "\nESTADISTICA DE SIMILITUD "
    "FINALIZADA CORRECTAMENTE"
)