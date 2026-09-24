import itertools
import numpy as np
import pandas as pd

# ============================================================
# DISPERSION MULTIVARIADA - DATASET PRIMARIO
# LINCS L1000
# 24 h - 10 uM
#
# Objetivo:
# comprobar si estradiol, calcitriol y calcipotriol
# presentan diferentes dispersiones multivariadas.
#
# Se calcula para cada firma la distancia Euclidea
# al centroide de su tratamiento.
#
# La significacion se evalua mediante permutacion EXACTA
# de las etiquetas de tratamiento dentro de cada bloque.
#
# 4 bloques x 3 tratamientos
# 3!^4 = 1296 permutaciones exactas
# ============================================================

ARCHIVO_EXPRESION = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_METADATA = "dataset_primario_metadata.csv"
ARCHIVO_GENES = "dataset_primario_genes.csv"

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol",
]

print("=" * 72)
print("DISPERSION MULTIVARIADA - DATASET PRIMARIO")
print("=" * 72)

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------

npz = np.load(ARCHIVO_EXPRESION)

if "expresion" in npz.files:
    X = npz["expresion"]

elif "expression" in npz.files:
    X = npz["expression"]

else:
    raise ValueError(
        "No se encontro la matriz de expresion."
    )

meta = pd.read_csv(ARCHIVO_METADATA)
genes = pd.read_csv(ARCHIVO_GENES)

print("\nForma de la matriz:")
print(X.shape)

# ------------------------------------------------------------
# 2. Controles basicos
# ------------------------------------------------------------

if X.shape[0] != 12:
    raise ValueError(
        "Se esperaban exactamente 12 firmas."
    )

if len(meta) != 12:
    raise ValueError(
        "La metadata no contiene exactamente 12 firmas."
    )

if X.shape[1] != len(genes):
    raise ValueError(
        "El numero de genes no coincide con la metadata."
    )

if "block_id" in meta.columns:
    BLOQUE = "block_id"
else:
    BLOQUE = "cell_id"

print("\nVariable utilizada como bloque:")
print(BLOQUE)

# ------------------------------------------------------------
# 3. Verificar diseño 4 x 3
# ------------------------------------------------------------

tabla = pd.crosstab(
    meta[BLOQUE],
    meta["pert_iname"]
)

print("\nDiseño bloque x tratamiento:")
print(tabla.to_string())

if tabla.shape != (4, 3):
    raise ValueError(
        "El diseño no contiene 4 bloques x 3 tratamientos."
    )

if not (tabla == 1).all().all():
    raise ValueError(
        "El diseño no contiene exactamente "
        "una firma por tratamiento y bloque."
    )

print("\nDiseño pareado: OK")

# ------------------------------------------------------------
# 4. Distancia de cada firma al centroide
# ------------------------------------------------------------

def distancias_a_centroides(Xdatos, etiquetas):
    """
    Calcula la distancia Euclidea de cada observacion
    al centroide del tratamiento al cual pertenece.
    """

    etiquetas = np.asarray(etiquetas)

    distancias = np.zeros(
        Xdatos.shape[0],
        dtype=float
    )

    for tratamiento in TRATAMIENTOS:

        idx = np.where(
            etiquetas == tratamiento
        )[0]

        if len(idx) == 0:
            raise ValueError(
                f"No hay observaciones para {tratamiento}"
            )

        centroide = np.mean(
            Xdatos[idx, :],
            axis=0
        )

        diferencias = (
            Xdatos[idx, :]
            - centroide
        )

        distancias[idx] = np.sqrt(
            np.sum(
                diferencias ** 2,
                axis=1
            )
        )

    return distancias


# ------------------------------------------------------------
# 5. Estadistica F de dispersion
# ------------------------------------------------------------

def estadistica_F_dispersion(
    Xdatos,
    etiquetas
):
    """
    ANOVA de una via aplicado a las distancias
    al centroide de cada tratamiento.

    La significacion NO se obtiene de la
    distribucion F teorica sino mediante
    permutacion exacta restringida por bloque.
    """

    etiquetas = np.asarray(etiquetas)

    distancias = distancias_a_centroides(
        Xdatos,
        etiquetas
    )

    media_global = np.mean(
        distancias
    )

    ss_between = 0.0
    ss_within = 0.0

    for tratamiento in TRATAMIENTOS:

        idx = np.where(
            etiquetas == tratamiento
        )[0]

        valores = distancias[idx]

        media_grupo = np.mean(
            valores
        )

        ss_between += (
            len(valores)
            * (media_grupo - media_global) ** 2
        )

        ss_within += np.sum(
            (valores - media_grupo) ** 2
        )

    k = len(TRATAMIENTOS)
    n = len(etiquetas)

    df_between = k - 1
    df_within = n - k

    ms_between = (
        ss_between / df_between
    )

    ms_within = (
        ss_within / df_within
    )

    if ms_within == 0:
        F = np.inf
    else:
        F = (
            ms_between / ms_within
        )

    return {
        "F": float(F),
        "df_between": int(df_between),
        "df_within": int(df_within),
        "distancias": distancias,
        "mean_distance":
            float(np.mean(distancias)),
    }


# ------------------------------------------------------------
# 6. Generar las 1296 permutaciones exactas
# ------------------------------------------------------------

def generar_permutaciones_exactas(
    meta,
    bloque_col
):

    bloques = list(
        pd.unique(meta[bloque_col])
    )

    indices_por_bloque = []
    permutaciones_por_bloque = []

    for bloque in bloques:

        indices = np.where(
            meta[bloque_col]
            .astype(str)
            .to_numpy()
            == str(bloque)
        )[0]

        etiquetas = (
            meta.iloc[indices]["pert_iname"]
            .astype(str)
            .tolist()
        )

        if len(indices) != 3:
            raise ValueError(
                f"El bloque {bloque} no tiene "
                "exactamente 3 firmas."
            )

        indices_por_bloque.append(
            indices
        )

        permutaciones_por_bloque.append(
            list(
                itertools.permutations(
                    etiquetas
                )
            )
        )

    for combinacion in itertools.product(
        *permutaciones_por_bloque
    ):

        perm = (
            meta["pert_iname"]
            .astype(str)
            .to_numpy()
            .copy()
        )

        for indices, etiquetas in zip(
            indices_por_bloque,
            combinacion
        ):

            perm[indices] = etiquetas

        yield perm


# ------------------------------------------------------------
# 7. Analisis completo
# ------------------------------------------------------------

def ejecutar_dispersion(
    Xdatos,
    nombre
):

    print("\n" + "=" * 72)
    print(nombre)
    print("=" * 72)

    etiquetas_obs = (
        meta["pert_iname"]
        .astype(str)
        .to_numpy()
    )

    observado = estadistica_F_dispersion(
        Xdatos,
        etiquetas_obs
    )

    F_obs = observado["F"]

    print("\nF observado:")
    print(f"{F_obs:.6f}")

    # ----------------------------------------
    # Distancias observadas
    # ----------------------------------------

    tabla_distancias = meta[
        [
            "sig_id",
            "pert_iname",
            "cell_id",
        ]
    ].copy()

    tabla_distancias[
        "distance_to_centroid"
    ] = observado["distancias"]

    archivo_distancias = (
        "DISPERSION_distancias_"
        + nombre
        + ".csv"
    )

    tabla_distancias.to_csv(
        archivo_distancias,
        index=False
    )

    # ----------------------------------------
    # Resumen por tratamiento
    # ----------------------------------------

    resumen = (
        tabla_distancias
        .groupby("pert_iname")
        ["distance_to_centroid"]
        .agg(
            [
                "count",
                "mean",
                "std",
                "median",
                "min",
                "max",
            ]
        )
        .reset_index()
    )

    print("\nDistancia al centroide por tratamiento:")
    print(
        resumen.to_string(
            index=False
        )
    )

    archivo_resumen = (
        "DISPERSION_resumen_"
        + nombre
        + ".csv"
    )

    resumen.to_csv(
        archivo_resumen,
        index=False
    )

    # ----------------------------------------
    # Permutaciones exactas
    # ----------------------------------------

    F_perm = []

    contador = 0

    for etiquetas_perm in (
        generar_permutaciones_exactas(
            meta,
            BLOQUE
        )
    ):

        resultado_perm = (
            estadistica_F_dispersion(
                Xdatos,
                etiquetas_perm
            )
        )

        F_perm.append(
            resultado_perm["F"]
        )

        contador += 1

    F_perm = np.asarray(
        F_perm,
        dtype=float
    )

    print("\nNumero de permutaciones exactas:")
    print(contador)

    if contador != 1296:
        raise ValueError(
            f"Se esperaban 1296 permutaciones "
            f"y se obtuvieron {contador}"
        )

    # La identidad esta incluida
    p_exacto = np.mean(
        F_perm >= (F_obs - 1e-12)
    )

    print("\np exacto:")
    print(f"{p_exacto:.8f}")

    # Guardar distribucion nula
    pd.DataFrame(
        {
            "F_permuted": F_perm
        }
    ).to_csv(
        "DISPERSION_permutaciones_"
        + nombre
        + ".csv",
        index=False
    )

    return {
        "analysis": nombre,
        "n_signatures": Xdatos.shape[0],
        "n_features": Xdatos.shape[1],
        "F_dispersion": F_obs,
        "df_between":
            observado["df_between"],
        "df_within":
            observado["df_within"],
        "n_exact_permutations":
            contador,
        "p_exact_dispersion":
            p_exacto,
    }


# ------------------------------------------------------------
# 8. Espacio completo: 12.328 genes
# ------------------------------------------------------------

resultado_full = ejecutar_dispersion(
    X,
    "full_12328_genes"
)

# ------------------------------------------------------------
# 9. Landmark: 978 genes
# ------------------------------------------------------------

genes["pr_is_lm"] = pd.to_numeric(
    genes["pr_is_lm"],
    errors="coerce"
)

mask_landmark = (
    genes["pr_is_lm"]
    .to_numpy()
    == 1
)

X_landmark = X[
    :,
    mask_landmark
]

print("\nGenes landmark:")
print(X_landmark.shape[1])

if X_landmark.shape[1] != 978:
    raise ValueError(
        "Se esperaban exactamente "
        "978 genes landmark."
    )

resultado_landmark = ejecutar_dispersion(
    X_landmark,
    "landmark_978_genes"
)

# ------------------------------------------------------------
# 10. Tabla final
# ------------------------------------------------------------

resultados = pd.DataFrame(
    [
        resultado_full,
        resultado_landmark,
    ]
)

resultados.to_csv(
    "DISPERSION_resultados_dataset_primario.csv",
    index=False
)

# ------------------------------------------------------------
# 11. Informe metodologico
# ------------------------------------------------------------

with open(
    "DISPERSION_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "ANALISIS DE DISPERSION MULTIVARIADA\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "Dataset primario:\n"
    )

    f.write(
        "12 firmas; 4 bloques completos; "
        "3 tratamientos por bloque.\n"
    )

    f.write(
        "Tiempo: 24 h.\n"
    )

    f.write(
        "Dosis: 10 uM.\n\n"
    )

    f.write(
        "Para cada tratamiento se calculo "
        "la distancia Euclidea de cada firma "
        "a su centroide multivariado.\n"
    )

    f.write(
        "Las diferencias entre tratamientos "
        "se cuantificaron mediante una "
        "estadistica F aplicada a dichas "
        "distancias.\n"
    )

    f.write(
        "La significacion se obtuvo mediante "
        "enumeracion exhaustiva de las "
        "1296 permutaciones posibles de las "
        "etiquetas de tratamiento dentro "
        "de los cuatro bloques experimentales.\n\n"
    )

    f.write(
        "El analisis se realizo tanto en el "
        "espacio completo de 12,328 genes "
        "como utilizando exclusivamente "
        "los 978 genes landmark.\n"
    )

print("\n" + "=" * 72)

print("Archivos creados:")

print(
    " - DISPERSION_resultados_dataset_primario.csv"
)

print(
    " - DISPERSION_resumen_full_12328_genes.csv"
)

print(
    " - DISPERSION_resumen_landmark_978_genes.csv"
)

print(
    " - DISPERSION_distancias_full_12328_genes.csv"
)

print(
    " - DISPERSION_distancias_landmark_978_genes.csv"
)

print(
    " - DISPERSION_permutaciones_full_12328_genes.csv"
)

print(
    " - DISPERSION_permutaciones_landmark_978_genes.csv"
)

print(
    " - DISPERSION_metodologia.txt"
)

print(
    "\nDISPERSION MULTIVARIADA FINALIZADA CORRECTAMENTE"
)