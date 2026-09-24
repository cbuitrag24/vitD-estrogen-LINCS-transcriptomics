import numpy as np
import pandas as pd

# ============================================================
# PROGRAMAS GENICOS CONSISTENTES
# Dataset primario LINCS L1000
# 24 h - 10 uM
#
# Objetivo:
# identificar genes cuya direccion de respuesta sea
# consistente entre los cuatro contextos celulares.
#
# No se aplican pruebas gen-a-gen con n=4.
# Se utiliza consistencia direccional + magnitud MODZ.
# ============================================================

ARCHIVO_EXPRESION = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_METADATA = "dataset_primario_metadata.csv"
ARCHIVO_GENES = "dataset_primario_genes.csv"

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

print("=" * 76)
print("PROGRAMAS GENICOS CONSISTENTES")
print("=" * 76)

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
        "No se encontro la matriz de expresion."
    )

meta = pd.read_csv(ARCHIVO_METADATA)
genes = pd.read_csv(ARCHIVO_GENES)

if X.shape != (12, 12328):
    raise ValueError(
        f"Forma inesperada: {X.shape}"
    )

if len(meta) != 12:
    raise ValueError(
        "Se esperaban 12 firmas."
    )

if len(genes) != 12328:
    raise ValueError(
        "Se esperaban 12328 genes."
    )

print("\nMatriz:")
print(X.shape)

# ------------------------------------------------------------
# 2. Preparar información génica
# ------------------------------------------------------------

genes["pr_is_lm"] = pd.to_numeric(
    genes["pr_is_lm"],
    errors="coerce"
).fillna(0).astype(int)

if genes["pr_is_lm"].sum() != 978:
    raise ValueError(
        "No se detectaron exactamente 978 genes landmark."
    )

# ------------------------------------------------------------
# 3. Función para resumir un tratamiento
# ------------------------------------------------------------

def resumir_tratamiento(tratamiento):

    idx = np.where(
        meta["pert_iname"]
        .astype(str)
        .to_numpy()
        == tratamiento
    )[0]

    if len(idx) != 4:
        raise ValueError(
            f"{tratamiento}: se esperaban 4 firmas "
            f"y se encontraron {len(idx)}."
        )

    Xt = X[idx, :]

    # Número de contextos con MODZ positivo/negativo
    n_up = np.sum(Xt > 0, axis=0)
    n_down = np.sum(Xt < 0, axis=0)
    n_zero = np.sum(Xt == 0, axis=0)

    mean_modz = np.mean(Xt, axis=0)
    median_modz = np.median(Xt, axis=0)

    mean_abs_modz = np.mean(
        np.abs(Xt),
        axis=0
    )

    # --------------------------------------------------------
    # Clasificación de consistencia
    # --------------------------------------------------------

    direction = np.full(
        X.shape[1],
        "context_dependent",
        dtype=object
    )

    consistency = np.full(
        X.shape[1],
        "2_or_less_of_4",
        dtype=object
    )

    # 3/4
    mask_up3 = n_up == 3
    mask_down3 = n_down == 3

    direction[mask_up3] = "up"
    direction[mask_down3] = "down"

    consistency[
        mask_up3 | mask_down3
    ] = "3_of_4"

    # 4/4
    mask_up4 = n_up == 4
    mask_down4 = n_down == 4

    direction[mask_up4] = "up"
    direction[mask_down4] = "down"

    consistency[
        mask_up4 | mask_down4
    ] = "4_of_4"

    salida = genes.copy()

    salida["treatment"] = tratamiento
    salida["n_up"] = n_up
    salida["n_down"] = n_down
    salida["n_zero"] = n_zero

    salida["direction"] = direction
    salida["consistency"] = consistency

    salida["mean_modz"] = mean_modz
    salida["median_modz"] = median_modz
    salida["mean_abs_modz"] = mean_abs_modz

    # Ranking:
    # primero consistencia y luego magnitud
    salida["consistency_rank"] = 0

    salida.loc[
        salida["consistency"] == "3_of_4",
        "consistency_rank"
    ] = 1

    salida.loc[
        salida["consistency"] == "4_of_4",
        "consistency_rank"
    ] = 2

    salida = salida.sort_values(
        [
            "consistency_rank",
            "mean_abs_modz"
        ],
        ascending=[
            False,
            False
        ]
    )

    return salida


# ------------------------------------------------------------
# 4. Ejecutar para los tres tratamientos
# ------------------------------------------------------------

todos = []

for tratamiento in TRATAMIENTOS:

    print("\n" + "-" * 76)
    print(tratamiento.upper())
    print("-" * 76)

    tabla = resumir_tratamiento(
        tratamiento
    )

    todos.append(tabla)

    archivo = (
        "PROGRAMAS_genes_"
        + tratamiento
        + ".csv"
    )

    tabla.to_csv(
        archivo,
        index=False
    )

    # --------------------------------------------------------
    # Resumen
    # --------------------------------------------------------

    resumen = (
        tabla
        .groupby(
            [
                "consistency",
                "direction"
            ]
        )
        .size()
        .reset_index(
            name="n_genes"
        )
    )

    print("\nResumen:")
    print(
        resumen.to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # Top genes 4/4
    # --------------------------------------------------------

    top = tabla[
        tabla["consistency"]
        == "4_of_4"
    ].copy()

    print("\nGenes 4/4 más fuertes:")

    columnas_top = [
        "pr_gene_symbol",
        "direction",
        "mean_modz",
        "median_modz",
        "mean_abs_modz",
        "pr_is_lm"
    ]

    print(
        top[
            columnas_top
        ]
        .head(20)
        .to_string(
            index=False
        )
    )

# ------------------------------------------------------------
# 5. Tabla conjunta
# ------------------------------------------------------------

todos = pd.concat(
    todos,
    ignore_index=True
)

todos.to_csv(
    "PROGRAMAS_genicos_todos_tratamientos.csv",
    index=False
)

# ------------------------------------------------------------
# 6. Resumen cuantitativo
# ------------------------------------------------------------

resumen_general = (
    todos
    .groupby(
        [
            "treatment",
            "consistency",
            "direction"
        ]
    )
    .size()
    .reset_index(
        name="n_genes"
    )
)

resumen_general.to_csv(
    "PROGRAMAS_resumen_consistencia.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Consenso 3/4 o 4/4
# ------------------------------------------------------------

consenso = todos[
    todos["consistency"]
    .isin(
        [
            "3_of_4",
            "4_of_4"
        ]
    )
].copy()

consenso.to_csv(
    "PROGRAMAS_genes_consenso_3o4.csv",
    index=False
)

# ------------------------------------------------------------
# 8. Consenso estricto 4/4
# ------------------------------------------------------------

consenso_estricto = todos[
    todos["consistency"]
    == "4_of_4"
].copy()

consenso_estricto.to_csv(
    "PROGRAMAS_genes_consenso_4o4.csv",
    index=False
)

# ------------------------------------------------------------
# 9. Resumen landmark
# ------------------------------------------------------------

landmark = todos[
    todos["pr_is_lm"]
    == 1
].copy()

resumen_landmark = (
    landmark
    .groupby(
        [
            "treatment",
            "consistency",
            "direction"
        ]
    )
    .size()
    .reset_index(
        name="n_genes"
    )
)

resumen_landmark.to_csv(
    "PROGRAMAS_resumen_landmark.csv",
    index=False
)

# ------------------------------------------------------------
# 10. Tabla compacta para inspección
# ------------------------------------------------------------

print("\n")
print("=" * 76)
print("RESUMEN GENERAL")
print("=" * 76)

print(
    resumen_general.to_string(
        index=False
    )
)

print("\n")
print("=" * 76)
print("RESUMEN LANDMARK")
print("=" * 76)

print(
    resumen_landmark.to_string(
        index=False
    )
)

# ------------------------------------------------------------
# 11. Informe metodológico
# ------------------------------------------------------------

with open(
    "PROGRAMAS_genicos_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "PROGRAMAS GENICOS CONSISTENTES\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "Se analizaron cuatro firmas "
        "independientes correspondientes "
        "a cuatro contextos celulares para "
        "cada tratamiento.\n\n"
    )

    f.write(
        "Para cada gen se contabilizó el "
        "número de contextos con MODZ positivo "
        "o negativo.\n\n"
    )

    f.write(
        "Los genes con dirección concordante "
        "en los cuatro contextos se "
        "clasificaron como 4_of_4.\n"
    )

    f.write(
        "Los genes concordantes en tres de "
        "cuatro contextos se clasificaron "
        "como 3_of_4.\n"
    )

    f.write(
        "Los restantes genes se consideraron "
        "dependientes del contexto para este "
        "análisis descriptivo.\n\n"
    )

    f.write(
        "La magnitud de respuesta se resumió "
        "mediante MODZ medio, MODZ mediano y "
        "MODZ absoluto medio.\n\n"
    )

    f.write(
        "No se utilizaron pruebas estadísticas "
        "gen-a-gen para definir los programas "
        "consenso debido al reducido número "
        "de contextos independientes.\n\n"
    )

    f.write(
        "Los genes landmark directamente "
        "medidos por L1000 fueron identificados "
        "para análisis de robustez.\n"
    )

print("\nArchivos creados:")
print(" - PROGRAMAS_genes_estradiol.csv")
print(" - PROGRAMAS_genes_calcitriol.csv")
print(" - PROGRAMAS_genes_calcipotriol.csv")
print(" - PROGRAMAS_genicos_todos_tratamientos.csv")
print(" - PROGRAMAS_resumen_consistencia.csv")
print(" - PROGRAMAS_genes_consenso_3o4.csv")
print(" - PROGRAMAS_genes_consenso_4o4.csv")
print(" - PROGRAMAS_resumen_landmark.csv")
print(" - PROGRAMAS_genicos_metodologia.txt")

print(
    "\nPROGRAMAS GENICOS "
    "FINALIZADOS CORRECTAMENTE"
)