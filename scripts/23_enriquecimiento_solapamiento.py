import numpy as np
import pandas as pd
from scipy.stats import hypergeom

# ============================================================
# ENRIQUECIMIENTO DE SOLAPAMIENTO ENTRE PROGRAMAS GENICOS
#
# Dataset primario LINCS L1000
# 24 h - 10 uM
#
# Programas:
# consenso estricto 4/4
#
# Se analiza por separado:
#   up
#   down
#
# Metricas:
#   interseccion observada
#   interseccion esperada por azar
#   Jaccard
#   fold enrichment
#   p hipergeometrico unilateral
#   FDR Benjamini-Hochberg
#
# Universos:
#   12,328 genes
#   978 genes landmark
# ============================================================

ARCHIVO = "PROGRAMAS_genicos_todos_tratamientos.csv"

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

PARES = [
    ("calcitriol", "calcipotriol"),
    ("estradiol", "calcitriol"),
    ("estradiol", "calcipotriol")
]

DIRECCIONES = [
    "up",
    "down"
]

print("=" * 80)
print("ENRIQUECIMIENTO DE SOLAPAMIENTO")
print("=" * 80)

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------

datos = pd.read_csv(ARCHIVO)

columnas_necesarias = [
    "pr_gene_id",
    "pr_gene_symbol",
    "pr_is_lm",
    "treatment",
    "direction",
    "consistency"
]

faltantes = [
    c for c in columnas_necesarias
    if c not in datos.columns
]

if faltantes:
    raise ValueError(
        f"Faltan columnas necesarias: {faltantes}"
    )

# ------------------------------------------------------------
# 2. Validar universo
# ------------------------------------------------------------

universo_total = (
    datos[
        ["pr_gene_id", "pr_is_lm"]
    ]
    .drop_duplicates(
        subset=["pr_gene_id"]
    )
    .copy()
)

N_FULL = len(universo_total)

N_LM = int(
    pd.to_numeric(
        universo_total["pr_is_lm"],
        errors="coerce"
    )
    .fillna(0)
    .sum()
)

print("\nUniverso completo:")
print(N_FULL)

print("\nUniverso landmark:")
print(N_LM)

if N_FULL != 12328:
    raise ValueError(
        f"Se esperaban 12328 genes y se encontraron {N_FULL}"
    )

if N_LM != 978:
    raise ValueError(
        f"Se esperaban 978 landmark y se encontraron {N_LM}"
    )

# ------------------------------------------------------------
# 3. Funcion BH
# ------------------------------------------------------------

def benjamini_hochberg(pvalues):

    pvalues = np.asarray(
        pvalues,
        dtype=float
    )

    n = len(pvalues)

    orden = np.argsort(
        pvalues
    )

    p_ordenados = (
        pvalues[orden]
    )

    q = np.empty(
        n,
        dtype=float
    )

    q_ordenados = (
        p_ordenados
        * n
        / np.arange(
            1,
            n + 1
        )
    )

    # Monotonicidad desde el final
    q_ordenados = np.minimum.accumulate(
        q_ordenados[::-1]
    )[::-1]

    q_ordenados = np.minimum(
        q_ordenados,
        1.0
    )

    q[orden] = q_ordenados

    return q

# ------------------------------------------------------------
# 4. Crear conjuntos
# ------------------------------------------------------------

def obtener_set(
    tratamiento,
    direccion,
    landmark_only=False
):

    sub = datos[
        (
            datos["treatment"]
            == tratamiento
        )
        &
        (
            datos["consistency"]
            == "4_of_4"
        )
        &
        (
            datos["direction"]
            == direccion
        )
    ].copy()

    if landmark_only:

        sub = sub[
            pd.to_numeric(
                sub["pr_is_lm"],
                errors="coerce"
            ).fillna(0)
            == 1
        ]

    return set(
        sub["pr_gene_id"]
        .astype(str)
    )

# ------------------------------------------------------------
# 5. Analisis de un par
# ------------------------------------------------------------

def analizar_par(
    tratamiento1,
    tratamiento2,
    direccion,
    universo,
    landmark_only,
    nombre_espacio
):

    A = obtener_set(
        tratamiento1,
        direccion,
        landmark_only
    )

    B = obtener_set(
        tratamiento2,
        direccion,
        landmark_only
    )

    n_A = len(A)
    n_B = len(B)

    interseccion = (
        A.intersection(B)
    )

    union = (
        A.union(B)
    )

    n_inter = len(
        interseccion
    )

    n_union = len(
        union
    )

    esperado = (
        n_A
        * n_B
        / universo
    )

    if esperado > 0:
        fold = (
            n_inter
            / esperado
        )
    else:
        fold = np.nan

    if n_union > 0:
        jaccard = (
            n_inter
            / n_union
        )
    else:
        jaccard = np.nan

    # --------------------------------------------------------
    # Hipergeometrica
    #
    # Universo = M
    # exitos en universo = n_A
    # extracciones = n_B
    # observado = n_inter
    #
    # P(X >= observado)
    # --------------------------------------------------------

    p_enrichment = hypergeom.sf(
        n_inter - 1,
        universo,
        n_A,
        n_B
    )

    return {
        "analysis_space":
            nombre_espacio,

        "comparison":
            tratamiento1
            + "_vs_"
            + tratamiento2,

        "treatment_1":
            tratamiento1,

        "treatment_2":
            tratamiento2,

        "direction":
            direccion,

        "universe_size":
            universo,

        "n_treatment_1":
            n_A,

        "n_treatment_2":
            n_B,

        "observed_overlap":
            n_inter,

        "expected_overlap":
            esperado,

        "fold_enrichment":
            fold,

        "jaccard":
            jaccard,

        "p_hypergeom":
            p_enrichment
    }

# ------------------------------------------------------------
# 6. Ejecutar todos los contrastes
# ------------------------------------------------------------

resultados = []

for t1, t2 in PARES:

    for direccion in DIRECCIONES:

        # Completo
        resultados.append(
            analizar_par(
                t1,
                t2,
                direccion,
                N_FULL,
                False,
                "full_12328_genes"
            )
        )

        # Landmark
        resultados.append(
            analizar_par(
                t1,
                t2,
                direccion,
                N_LM,
                True,
                "landmark_978_genes"
            )
        )

resultados = pd.DataFrame(
    resultados
)

# ------------------------------------------------------------
# 7. Correccion FDR
#
# Se corrigen conjuntamente los 12 contrastes.
# ------------------------------------------------------------

resultados[
    "q_BH"
] = benjamini_hochberg(
    resultados["p_hypergeom"]
    .to_numpy()
)

# ------------------------------------------------------------
# 8. Guardar tabla
# ------------------------------------------------------------

resultados.to_csv(
    "ENRIQUECIMIENTO_solapamiento_resultados.csv",
    index=False
)

# ------------------------------------------------------------
# 9. Resumen completo
# ------------------------------------------------------------

print("\n")
print("=" * 80)
print("RESULTADOS")
print("=" * 80)

columnas_mostrar = [
    "analysis_space",
    "comparison",
    "direction",
    "n_treatment_1",
    "n_treatment_2",
    "observed_overlap",
    "expected_overlap",
    "fold_enrichment",
    "jaccard",
    "p_hypergeom",
    "q_BH"
]

print(
    resultados[
        columnas_mostrar
    ]
    .round(
        {
            "expected_overlap": 3,
            "fold_enrichment": 3,
            "jaccard": 4,
            "p_hypergeom": 6,
            "q_BH": 6
        }
    )
    .to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 10. Resumen agregado same-direction
#
# Para comparar con los resultados previos:
# suma de solapamientos up + down.
# El p no se combina aquí.
# ------------------------------------------------------------

agregado = (
    resultados
    .groupby(
        [
            "analysis_space",
            "comparison"
        ]
    )
    .agg(
        observed_same_direction=(
            "observed_overlap",
            "sum"
        ),
        expected_same_direction=(
            "expected_overlap",
            "sum"
        )
    )
    .reset_index()
)

agregado[
    "fold_same_direction"
] = (
    agregado[
        "observed_same_direction"
    ]
    /
    agregado[
        "expected_same_direction"
    ]
)

agregado.to_csv(
    "ENRIQUECIMIENTO_resumen_same_direction.csv",
    index=False
)

print("\n")
print("=" * 80)
print("RESUMEN MISMA DIRECCION: UP + DOWN")
print("=" * 80)

print(
    agregado.round(
        {
            "expected_same_direction": 3,
            "fold_same_direction": 3
        }
    ).to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 11. Guardar genes compartidos
# ------------------------------------------------------------

for espacio, lm in [
    ("full_12328_genes", False),
    ("landmark_978_genes", True)
]:

    for t1, t2 in PARES:

        filas = []

        for direccion in DIRECCIONES:

            A = obtener_set(
                t1,
                direccion,
                lm
            )

            B = obtener_set(
                t2,
                direccion,
                lm
            )

            compartidos = (
                A.intersection(B)
            )

            for gene_id in compartidos:

                filas.append(
                    {
                        "analysis_space":
                            espacio,
                        "comparison":
                            t1
                            + "_vs_"
                            + t2,
                        "direction":
                            direccion,
                        "pr_gene_id":
                            gene_id
                    }
                )

        tabla_genes = pd.DataFrame(
            filas
        )

        if len(tabla_genes) > 0:

            info = (
                datos[
                    [
                        "pr_gene_id",
                        "pr_gene_symbol",
                        "pr_is_lm"
                    ]
                ]
                .drop_duplicates(
                    subset=[
                        "pr_gene_id"
                    ]
                )
                .copy()
            )

            info[
                "pr_gene_id"
            ] = info[
                "pr_gene_id"
            ].astype(str)

            tabla_genes = tabla_genes.merge(
                info,
                on="pr_gene_id",
                how="left"
            )

        tabla_genes.to_csv(
            "ENRIQUECIMIENTO_genes_compartidos_"
            + espacio
            + "_"
            + t1
            + "_vs_"
            + t2
            + ".csv",
            index=False
        )

# ------------------------------------------------------------
# 12. Informe metodologico
# ------------------------------------------------------------

with open(
    "ENRIQUECIMIENTO_solapamiento_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "ENRIQUECIMIENTO DE SOLAPAMIENTO "
        "ENTRE PROGRAMAS GENICOS\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "El análisis utilizó los programas "
        "direccionales estrictos 4/4.\n\n"
    )

    f.write(
        "Los genes regulados hacia arriba "
        "y hacia abajo se analizaron "
        "por separado.\n\n"
    )

    f.write(
        "Para cada par de tratamientos se "
        "calculó el número observado de "
        "genes compartidos en la misma "
        "dirección, el solapamiento esperado "
        "por azar, el índice de Jaccard y "
        "el enriquecimiento observado/esperado.\n\n"
    )

    f.write(
        "La significación del solapamiento "
        "se evaluó mediante una prueba "
        "hipergeométrica unilateral para "
        "enriquecimiento.\n\n"
    )

    f.write(
        "El universo primario estuvo formado "
        "por los 12,328 genes presentes en "
        "la matriz LINCS L1000 analizada.\n"
    )

    f.write(
        "Como análisis de sensibilidad, "
        "el procedimiento se repitió "
        "restringiendo el universo a los "
        "978 genes landmark directamente "
        "medidos por L1000.\n\n"
    )

    f.write(
        "Los valores p de los 12 contrastes "
        "de enriquecimiento se corrigieron "
        "mediante el procedimiento de "
        "Benjamini-Hochberg.\n"
    )

print("\nArchivos creados:")
print(
    " - ENRIQUECIMIENTO_solapamiento_resultados.csv"
)
print(
    " - ENRIQUECIMIENTO_resumen_same_direction.csv"
)
print(
    " - archivos ENRIQUECIMIENTO_genes_compartidos_*.csv"
)
print(
    " - ENRIQUECIMIENTO_solapamiento_metodologia.txt"
)

print(
    "\nENRIQUECIMIENTO DE SOLAPAMIENTO "
    "FINALIZADO CORRECTAMENTE"
)