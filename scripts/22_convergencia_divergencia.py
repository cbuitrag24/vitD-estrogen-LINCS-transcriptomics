import numpy as np
import pandas as pd

# ============================================================
# CONVERGENCIA Y DIVERGENCIA GENICA
# Dataset primario LINCS L1000
# 24 h - 10 uM
#
# Se utilizan exclusivamente programas direccionales
# estrictos 4/4 para definir:
#
# 1. Convergencia = mismo gen, misma direccion
# 2. Divergencia   = mismo gen, direccion opuesta
# 3. Nucleo comun a los tres tratamientos
# 4. Nucleo compartido calcitriol-calcipotriol
# 5. Estado del nucleo VitD frente a estradiol
# 6. Programas estrictos 4/4 exclusivos
#
# IMPORTANTE:
# "exclusivo" significa exclusivo entre los programas
# estrictos 4/4, no ausencia completa de respuesta.
# ============================================================

ARCHIVO = "PROGRAMAS_genicos_todos_tratamientos.csv"

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

print("=" * 78)
print("CONVERGENCIA Y DIVERGENCIA GENICA")
print("=" * 78)

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
    "consistency",
    "mean_modz",
    "median_modz",
    "mean_abs_modz"
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
# 2. Seleccionar consenso estricto 4/4
# ------------------------------------------------------------

estricto = datos[
    datos["consistency"] == "4_of_4"
].copy()

print("\nGenes estrictos 4/4 por tratamiento:")

print(
    estricto
    .groupby(
        ["treatment", "direction"]
    )
    .size()
    .reset_index(name="n_genes")
    .to_string(index=False)
)

# ------------------------------------------------------------
# 3. Crear una tabla de estado para todos los genes
# ------------------------------------------------------------

base_genes = (
    datos[
        [
            "pr_gene_id",
            "pr_gene_symbol",
            "pr_is_lm"
        ]
    ]
    .drop_duplicates(
        subset=["pr_gene_id"]
    )
    .copy()
)

# Para cada tratamiento:
# strict = 1/0
# direction
# mean_modz

for tratamiento in TRATAMIENTOS:

    sub = datos[
        datos["treatment"] == tratamiento
    ][
        [
            "pr_gene_id",
            "direction",
            "consistency",
            "mean_modz",
            "median_modz",
            "mean_abs_modz"
        ]
    ].copy()

    sub[
        f"{tratamiento}_strict4"
    ] = (
        sub["consistency"]
        == "4_of_4"
    ).astype(int)

    sub = sub.rename(
        columns={
            "direction":
                f"{tratamiento}_direction",
            "consistency":
                f"{tratamiento}_consistency",
            "mean_modz":
                f"{tratamiento}_mean_modz",
            "median_modz":
                f"{tratamiento}_median_modz",
            "mean_abs_modz":
                f"{tratamiento}_mean_abs_modz"
        }
    )

    columnas = [
        "pr_gene_id",
        f"{tratamiento}_direction",
        f"{tratamiento}_consistency",
        f"{tratamiento}_mean_modz",
        f"{tratamiento}_median_modz",
        f"{tratamiento}_mean_abs_modz",
        f"{tratamiento}_strict4"
    ]

    base_genes = base_genes.merge(
        sub[columnas],
        on="pr_gene_id",
        how="left"
    )

# ------------------------------------------------------------
# 4. Función para clasificar un par de tratamientos
# ------------------------------------------------------------

def analizar_par(
    tabla,
    tratamiento1,
    tratamiento2
):

    strict1 = (
        tabla[f"{tratamiento1}_strict4"]
        == 1
    )

    strict2 = (
        tabla[f"{tratamiento2}_strict4"]
        == 1
    )

    ambos = strict1 & strict2

    misma_direccion = (
        tabla[f"{tratamiento1}_direction"]
        ==
        tabla[f"{tratamiento2}_direction"]
    )

    convergente = (
        ambos
        & misma_direccion
    )

    divergente = (
        ambos
        & (~misma_direccion)
    )

    solo_t1 = (
        strict1
        & (~strict2)
    )

    solo_t2 = (
        strict2
        & (~strict1)
    )

    categoria = np.full(
        len(tabla),
        "neither_strict4",
        dtype=object
    )

    categoria[solo_t1] = (
        f"strict4_only_{tratamiento1}"
    )

    categoria[solo_t2] = (
        f"strict4_only_{tratamiento2}"
    )

    categoria[convergente] = (
        "convergent_same_direction"
    )

    categoria[divergente] = (
        "divergent_opposite_direction"
    )

    salida = tabla[
        [
            "pr_gene_id",
            "pr_gene_symbol",
            "pr_is_lm",
            f"{tratamiento1}_direction",
            f"{tratamiento2}_direction",
            f"{tratamiento1}_mean_modz",
            f"{tratamiento2}_mean_modz",
            f"{tratamiento1}_mean_abs_modz",
            f"{tratamiento2}_mean_abs_modz"
        ]
    ].copy()

    salida["treatment_1"] = tratamiento1
    salida["treatment_2"] = tratamiento2
    salida["category"] = categoria

    # Magnitud media conjunta para ranking descriptivo
    salida["joint_mean_abs_modz"] = (
        (
            salida[
                f"{tratamiento1}_mean_abs_modz"
            ]
            +
            salida[
                f"{tratamiento2}_mean_abs_modz"
            ]
        )
        / 2.0
    )

    return salida


# ------------------------------------------------------------
# 5. Analizar los tres pares
# ------------------------------------------------------------

pares = [
    ("calcitriol", "calcipotriol"),
    ("estradiol", "calcitriol"),
    ("estradiol", "calcipotriol")
]

resultados_pares = []
resumen_pares = []

for t1, t2 in pares:

    salida = analizar_par(
        base_genes,
        t1,
        t2
    )

    resultados_pares.append(
        salida
    )

    resumen = (
        salida
        .groupby("category")
        .size()
        .reset_index(
            name="n_genes"
        )
    )

    resumen["comparison"] = (
        t1 + "_vs_" + t2
    )

    resumen_pares.append(
        resumen
    )

    # Guardar solamente genes que participan
    # en algún programa estricto del par
    relevante = salida[
        salida["category"]
        != "neither_strict4"
    ].copy()

    relevante = relevante.sort_values(
        "joint_mean_abs_modz",
        ascending=False
    )

    relevante.to_csv(
        "CONVERGENCIA_"
        + t1
        + "_vs_"
        + t2
        + ".csv",
        index=False
    )

# ------------------------------------------------------------
# 6. Tabla conjunta y resumen
# ------------------------------------------------------------

resultados_pares = pd.concat(
    resultados_pares,
    ignore_index=True
)

resumen_pares = pd.concat(
    resumen_pares,
    ignore_index=True
)

resultados_pares.to_csv(
    "CONVERGENCIA_todos_los_pares.csv",
    index=False
)

resumen_pares.to_csv(
    "CONVERGENCIA_resumen_pares.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Núcleo común a los tres tratamientos
# ------------------------------------------------------------

todos_strict = (
    (base_genes["estradiol_strict4"] == 1)
    &
    (base_genes["calcitriol_strict4"] == 1)
    &
    (base_genes["calcipotriol_strict4"] == 1)
)

misma_direccion_tres = (
    (
        base_genes["estradiol_direction"]
        ==
        base_genes["calcitriol_direction"]
    )
    &
    (
        base_genes["estradiol_direction"]
        ==
        base_genes["calcipotriol_direction"]
    )
)

nucleo_tres_convergente = base_genes[
    todos_strict
    &
    misma_direccion_tres
].copy()

nucleo_tres_divergente = base_genes[
    todos_strict
    &
    (~misma_direccion_tres)
].copy()

# Ranking por magnitud media de los 3 tratamientos
for tabla in [
    nucleo_tres_convergente,
    nucleo_tres_divergente
]:

    tabla["joint_mean_abs_modz_3"] = (
        (
            tabla["estradiol_mean_abs_modz"]
            +
            tabla["calcitriol_mean_abs_modz"]
            +
            tabla["calcipotriol_mean_abs_modz"]
        )
        / 3.0
    )

nucleo_tres_convergente = (
    nucleo_tres_convergente
    .sort_values(
        "joint_mean_abs_modz_3",
        ascending=False
    )
)

nucleo_tres_divergente = (
    nucleo_tres_divergente
    .sort_values(
        "joint_mean_abs_modz_3",
        ascending=False
    )
)

nucleo_tres_convergente.to_csv(
    "CONVERGENCIA_nucleo_tres_misma_direccion.csv",
    index=False
)

nucleo_tres_divergente.to_csv(
    "CONVERGENCIA_nucleo_tres_direccion_opuesta.csv",
    index=False
)

# ------------------------------------------------------------
# 8. Núcleo estricto vitamina D
#    calcitriol + calcipotriol
#    ambos 4/4 y misma dirección
# ------------------------------------------------------------

vitd_strict = (
    (base_genes["calcitriol_strict4"] == 1)
    &
    (base_genes["calcipotriol_strict4"] == 1)
    &
    (
        base_genes["calcitriol_direction"]
        ==
        base_genes["calcipotriol_direction"]
    )
)

nucleo_vitd = base_genes[
    vitd_strict
].copy()

# ------------------------------------------------------------
# 9. Clasificar cómo responde estradiol
#    frente al núcleo VitD
# ------------------------------------------------------------

estado_e2 = []

for _, fila in nucleo_vitd.iterrows():

    if fila["estradiol_strict4"] != 1:

        estado = "E2_not_strict4"

    else:

        direccion_vitd = (
            fila["calcitriol_direction"]
        )

        direccion_e2 = (
            fila["estradiol_direction"]
        )

        if direccion_e2 == direccion_vitd:
            estado = "E2_same_direction"
        else:
            estado = "E2_opposite_direction"

    estado_e2.append(
        estado
    )

nucleo_vitd[
    "E2_status_relative_to_VitD"
] = estado_e2

nucleo_vitd[
    "vitd_joint_mean_abs_modz"
] = (
    (
        nucleo_vitd[
            "calcitriol_mean_abs_modz"
        ]
        +
        nucleo_vitd[
            "calcipotriol_mean_abs_modz"
        ]
    )
    / 2.0
)

nucleo_vitd = nucleo_vitd.sort_values(
    "vitd_joint_mean_abs_modz",
    ascending=False
)

nucleo_vitd.to_csv(
    "CONVERGENCIA_nucleo_VitD_y_estado_E2.csv",
    index=False
)

resumen_vitd_e2 = (
    nucleo_vitd
    .groupby(
        "E2_status_relative_to_VitD"
    )
    .size()
    .reset_index(
        name="n_genes"
    )
)

resumen_vitd_e2.to_csv(
    "CONVERGENCIA_resumen_nucleo_VitD_vs_E2.csv",
    index=False
)

# ------------------------------------------------------------
# 10. Programas estrictos exclusivos
#
# Exclusivo = 4/4 solamente en ese tratamiento
# respecto a los otros dos programas 4/4.
# ------------------------------------------------------------

exclusivos = []

for tratamiento in TRATAMIENTOS:

    otros = [
        x for x in TRATAMIENTOS
        if x != tratamiento
    ]

    mask = (
        (base_genes[
            f"{tratamiento}_strict4"
        ] == 1)
        &
        (base_genes[
            f"{otros[0]}_strict4"
        ] == 0)
        &
        (base_genes[
            f"{otros[1]}_strict4"
        ] == 0)
    )

    temp = base_genes[
        mask
    ].copy()

    temp["strict4_exclusive_to"] = (
        tratamiento
    )

    temp = temp.sort_values(
        f"{tratamiento}_mean_abs_modz",
        ascending=False
    )

    temp.to_csv(
        "CONVERGENCIA_programa_estricto_exclusivo_"
        + tratamiento
        + ".csv",
        index=False
    )

    exclusivos.append(
        temp
    )

# ------------------------------------------------------------
# 11. Resumen principal
# ------------------------------------------------------------

print("\n")
print("=" * 78)
print("RESUMEN DE COMPARACIONES PAREADAS")
print("=" * 78)

resumen_importante = resumen_pares[
    resumen_pares["category"].isin(
        [
            "convergent_same_direction",
            "divergent_opposite_direction"
        ]
    )
].copy()

print(
    resumen_importante[
        [
            "comparison",
            "category",
            "n_genes"
        ]
    ].to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 12. Núcleo común a los tres
# ------------------------------------------------------------

print("\n")
print("=" * 78)
print("NUCLEO ESTRICTO COMPARTIDO POR LOS TRES")
print("=" * 78)

print(
    "Misma direccion:",
    len(nucleo_tres_convergente)
)

print(
    "Direcciones opuestas/mixtas:",
    len(nucleo_tres_divergente)
)

# ------------------------------------------------------------
# 13. Núcleo vitamina D y relación con E2
# ------------------------------------------------------------

print("\n")
print("=" * 78)
print("NUCLEO ESTRICTO VITAMINA D")
print("=" * 78)

print(
    "Genes calcitriol-calcipotriol "
    "4/4 en la misma direccion:",
    len(nucleo_vitd)
)

print("\nEstado de esos genes frente a estradiol:")

print(
    resumen_vitd_e2.to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 14. Exclusivos
# ------------------------------------------------------------

print("\n")
print("=" * 78)
print("PROGRAMAS ESTRICTOS 4/4 EXCLUSIVOS")
print("=" * 78)

for tratamiento, tabla in zip(
    TRATAMIENTOS,
    exclusivos
):

    print(
        f"{tratamiento}: {len(tabla)} genes"
    )

# ------------------------------------------------------------
# 15. Landmark dentro del núcleo VitD
# ------------------------------------------------------------

nucleo_vitd_landmark = nucleo_vitd[
    nucleo_vitd["pr_is_lm"] == 1
].copy()

nucleo_vitd_landmark.to_csv(
    "CONVERGENCIA_nucleo_VitD_landmark.csv",
    index=False
)

print("\nGenes landmark en núcleo VitD:")
print(
    len(nucleo_vitd_landmark)
)

# ------------------------------------------------------------
# 16. Top 20 núcleo VitD
# ------------------------------------------------------------

print("\n")
print("=" * 78)
print("TOP 20 DEL NUCLEO VitD POR MAGNITUD")
print("=" * 78)

columnas_top = [
    "pr_gene_symbol",
    "calcitriol_direction",
    "calcitriol_mean_modz",
    "calcipotriol_mean_modz",
    "vitd_joint_mean_abs_modz",
    "E2_status_relative_to_VitD",
    "estradiol_mean_modz",
    "pr_is_lm"
]

print(
    nucleo_vitd[
        columnas_top
    ]
    .head(20)
    .round(4)
    .to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 17. Informe metodológico
# ------------------------------------------------------------

with open(
    "CONVERGENCIA_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "CONVERGENCIA Y DIVERGENCIA GENICA\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "El análisis se restringió inicialmente "
        "a genes con respuesta direccional "
        "concordante en los cuatro contextos "
        "celulares (consenso estricto 4/4).\n\n"
    )

    f.write(
        "Para cada par de tratamientos, "
        "la convergencia se definió como un "
        "gen perteneciente al consenso 4/4 "
        "de ambos tratamientos y regulado "
        "en la misma dirección.\n\n"
    )

    f.write(
        "La divergencia se definió como un "
        "gen perteneciente al consenso 4/4 "
        "de ambos tratamientos pero regulado "
        "en direcciones opuestas.\n\n"
    )

    f.write(
        "Se definió además un núcleo vitamina D "
        "formado por genes concordantes 4/4 "
        "tanto para calcitriol como para "
        "calcipotriol y regulados en la misma "
        "dirección por ambos compuestos.\n\n"
    )

    f.write(
        "Para cada gen del núcleo vitamina D "
        "se determinó si estradiol presentaba "
        "un consenso estricto 4/4 en la misma "
        "dirección, en dirección opuesta, "
        "o si no alcanzaba el criterio 4/4.\n\n"
    )

    f.write(
        "Los programas denominados exclusivos "
        "corresponden únicamente a exclusividad "
        "entre los consensos estrictos 4/4. "
        "La ausencia de un gen del consenso "
        "4/4 de otro tratamiento no implica "
        "ausencia de respuesta transcriptómica "
        "a dicho tratamiento.\n\n"
    )

    f.write(
        "Las magnitudes MODZ se utilizaron "
        "para ordenar descriptivamente los genes, "
        "no como criterio de significación "
        "estadística gen-a-gen.\n"
    )

print("\n")
print("=" * 78)
print("ARCHIVOS CREADOS")
print("=" * 78)

print(
    " - CONVERGENCIA_calcitriol_vs_calcipotriol.csv"
)
print(
    " - CONVERGENCIA_estradiol_vs_calcitriol.csv"
)
print(
    " - CONVERGENCIA_estradiol_vs_calcipotriol.csv"
)
print(
    " - CONVERGENCIA_todos_los_pares.csv"
)
print(
    " - CONVERGENCIA_resumen_pares.csv"
)
print(
    " - CONVERGENCIA_nucleo_tres_misma_direccion.csv"
)
print(
    " - CONVERGENCIA_nucleo_tres_direccion_opuesta.csv"
)
print(
    " - CONVERGENCIA_nucleo_VitD_y_estado_E2.csv"
)
print(
    " - CONVERGENCIA_resumen_nucleo_VitD_vs_E2.csv"
)
print(
    " - CONVERGENCIA_nucleo_VitD_landmark.csv"
)
print(
    " - CONVERGENCIA_programa_estricto_exclusivo_estradiol.csv"
)
print(
    " - CONVERGENCIA_programa_estricto_exclusivo_calcitriol.csv"
)
print(
    " - CONVERGENCIA_programa_estricto_exclusivo_calcipotriol.csv"
)
print(
    " - CONVERGENCIA_metodologia.txt"
)

print(
    "\nCONVERGENCIA Y DIVERGENCIA "
    "FINALIZADAS CORRECTAMENTE"
)