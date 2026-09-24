import itertools
import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform

# ============================================================
# PERMANOVA PAREADA - DATASET PRIMARIO
# LINCS L1000 - 24 h - 10 uM
#
# Diseño:
# 4 bloques experimentales completos
# 3 tratamientos por bloque:
#   estradiol
#   calcitriol
#   calcipotriol
#
# Prueba principal:
# efecto global del tratamiento
#
# Permutacion:
# etiquetas de tratamiento permutadas DENTRO de cada bloque
#
# Como existen 3! permutaciones por bloque y 4 bloques:
# 6^4 = 1296 permutaciones posibles.
#
# Se utiliza enumeracion EXACTA, no muestreo aleatorio.
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
print("PERMANOVA PAREADA EXACTA - DATASET PRIMARIO")
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
# 2. Controles
# ------------------------------------------------------------

if X.shape[0] != 12:
    raise ValueError(
        "Se esperaban exactamente 12 firmas."
    )

if len(meta) != 12:
    raise ValueError(
        "Metadata no contiene exactamente 12 firmas."
    )

if X.shape[1] != len(genes):
    raise ValueError(
        "Numero de genes no coincide con metadata."
    )

# Elegir identificador real del bloque
if "block_id" in meta.columns:

    BLOQUE = "block_id"

else:

    # En este dataset primario cada linea celular
    # corresponde al bloque experimental completo.
    BLOQUE = "cell_id"

print("\nVariable utilizada como bloque:")
print(BLOQUE)

# ------------------------------------------------------------
# 3. Verificar diseño completo 1:1:1
# ------------------------------------------------------------

tabla_diseno = pd.crosstab(
    meta[BLOQUE],
    meta["pert_iname"]
)

print("\nDiseño bloque x tratamiento:")
print(tabla_diseno.to_string())

if tabla_diseno.shape != (4, 3):
    raise ValueError(
        "El diseño no contiene 4 bloques x 3 tratamientos."
    )

if not (tabla_diseno == 1).all().all():
    raise ValueError(
        "El diseño no es exactamente 1:1:1."
    )

print("\nDiseño pareado: OK")

# ------------------------------------------------------------
# 4. Funciones para PERMANOVA
# ------------------------------------------------------------

def matriz_modelo(bloques, tratamientos):
    """
    Modelo aditivo:
    intercepto + bloque + tratamiento
    """

    bloques = np.asarray(bloques)
    tratamientos = np.asarray(tratamientos)

    categorias_bloque = sorted(
        pd.unique(bloques)
    )

    categorias_trat = TRATAMIENTOS

    columnas = [
        np.ones(len(bloques))
    ]

    # Dummy coding para bloque
    for categoria in categorias_bloque[1:]:

        columnas.append(
            (bloques == categoria).astype(float)
        )

    # Dummy coding para tratamiento
    for categoria in categorias_trat[1:]:

        columnas.append(
            (tratamientos == categoria).astype(float)
        )

    return np.column_stack(columnas)


def matriz_modelo_bloque(bloques):
    """
    Modelo reducido:
    intercepto + bloque
    """

    bloques = np.asarray(bloques)

    categorias = sorted(
        pd.unique(bloques)
    )

    columnas = [
        np.ones(len(bloques))
    ]

    for categoria in categorias[1:]:

        columnas.append(
            (bloques == categoria).astype(float)
        )

    return np.column_stack(columnas)


def matriz_modelo_tratamiento(tratamientos):
    """
    Modelo:
    intercepto + tratamiento
    """

    tratamientos = np.asarray(tratamientos)

    columnas = [
        np.ones(len(tratamientos))
    ]

    for categoria in TRATAMIENTOS[1:]:

        columnas.append(
            (tratamientos == categoria).astype(float)
        )

    return np.column_stack(columnas)


def hat_matrix(Xmodelo):
    """
    Matriz de proyeccion H.
    """

    return (
        Xmodelo
        @ np.linalg.pinv(Xmodelo)
    )


def gower_matrix(Xdatos):
    """
    Construye matriz de distancias Euclideas
    y luego matriz centrada de Gower.
    """

    D = squareform(
        pdist(
            Xdatos,
            metric="euclidean"
        )
    )

    n = D.shape[0]

    J = (
        np.eye(n)
        - np.ones((n, n)) / n
    )

    G = (
        -0.5
        * J
        @ (D ** 2)
        @ J
    )

    return D, G


def estadistica_tratamiento(
    G,
    bloques,
    tratamientos
):
    """
    Efecto marginal del tratamiento
    controlando por bloque.
    """

    X_full = matriz_modelo(
        bloques,
        tratamientos
    )

    X_red = matriz_modelo_bloque(
        bloques
    )

    H_full = hat_matrix(X_full)
    H_red = hat_matrix(X_red)

    n = len(tratamientos)

    rango_full = np.linalg.matrix_rank(
        X_full
    )

    rango_red = np.linalg.matrix_rank(
        X_red
    )

    df_efecto = (
        rango_full
        - rango_red
    )

    df_residual = (
        n
        - rango_full
    )

    SS_efecto = np.trace(
        (H_full - H_red) @ G
    )

    SS_residual = np.trace(
        (np.eye(n) - H_full) @ G
    )

    SS_total = np.trace(G)

    MS_efecto = (
        SS_efecto / df_efecto
    )

    MS_residual = (
        SS_residual / df_residual
    )

    F = (
        MS_efecto / MS_residual
    )

    R2_total = (
        SS_efecto / SS_total
    )

    R2_parcial = (
        SS_efecto
        / (SS_efecto + SS_residual)
    )

    return {
        "F": float(F),
        "SS_effect": float(SS_efecto),
        "SS_residual": float(SS_residual),
        "SS_total": float(SS_total),
        "df_effect": int(df_efecto),
        "df_residual": int(df_residual),
        "R2_total": float(R2_total),
        "partial_R2": float(R2_parcial),
    }


def varianza_bloque(
    G,
    bloques,
    tratamientos
):
    """
    Cuantifica la contribucion del bloque/contexto.

    IMPORTANTE:
    se informa descriptivamente.
    NO se asigna p inferencial porque en este
    dataset cada linea celular esta ligada a
    un unico bloque experimental.
    """

    X_full = matriz_modelo(
        bloques,
        tratamientos
    )

    X_red = matriz_modelo_tratamiento(
        tratamientos
    )

    H_full = hat_matrix(X_full)
    H_red = hat_matrix(X_red)

    SS_bloque = np.trace(
        (H_full - H_red) @ G
    )

    SS_total = np.trace(G)

    return {
        "SS_block_context":
            float(SS_bloque),

        "R2_block_context":
            float(SS_bloque / SS_total),
    }


# ------------------------------------------------------------
# 5. Permutaciones exactas dentro de bloque
# ------------------------------------------------------------

def generar_permutaciones_exactas(
    meta,
    bloque_col
):

    bloques = list(
        pd.unique(meta[bloque_col])
    )

    indices_por_bloque = []

    etiquetas_por_bloque = []

    for bloque in bloques:

        indices = np.where(
            meta[bloque_col].to_numpy()
            == bloque
        )[0]

        etiquetas = (
            meta.iloc[indices]["pert_iname"]
            .tolist()
        )

        if len(indices) != 3:
            raise ValueError(
                f"Bloque {bloque} no tiene 3 firmas."
            )

        indices_por_bloque.append(indices)

        etiquetas_por_bloque.append(
            list(
                itertools.permutations(
                    etiquetas
                )
            )
        )

    for combinacion in itertools.product(
        *etiquetas_por_bloque
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
# 6. Funcion completa para un espacio genetico
# ------------------------------------------------------------

def ejecutar_permanova(
    Xdatos,
    nombre_analisis
):

    print("\n" + "=" * 72)
    print(nombre_analisis)
    print("=" * 72)

    _, G = gower_matrix(
        Xdatos
    )

    bloques = (
        meta[BLOQUE]
        .astype(str)
        .to_numpy()
    )

    tratamientos_obs = (
        meta["pert_iname"]
        .astype(str)
        .to_numpy()
    )

    # Estadistica observada
    observado = estadistica_tratamiento(
        G,
        bloques,
        tratamientos_obs
    )

    contexto = varianza_bloque(
        G,
        bloques,
        tratamientos_obs
    )

    F_obs = observado["F"]

    print("\nPseudo-F observado para tratamiento:")
    print(f"{F_obs:.6f}")

    print("\nR2 tratamiento:")
    print(
        f"{observado['R2_total']:.6f}"
    )

    print("\nR2 parcial tratamiento:")
    print(
        f"{observado['partial_R2']:.6f}"
    )

    print("\nR2 bloque/contexto:")
    print(
        f"{contexto['R2_block_context']:.6f}"
    )

    # ----------------------------------------
    # Permutacion exacta
    # ----------------------------------------

    valores_F = []

    contador = 0

    for perm in generar_permutaciones_exactas(
        meta,
        BLOQUE
    ):

        resultado_perm = (
            estadistica_tratamiento(
                G,
                bloques,
                perm
            )
        )

        valores_F.append(
            resultado_perm["F"]
        )

        contador += 1

    valores_F = np.asarray(
        valores_F,
        dtype=float
    )

    print("\nNumero de permutaciones exactas:")
    print(contador)

    if contador != 1296:

        raise ValueError(
            f"Se esperaban 1296 permutaciones "
            f"y se obtuvieron {contador}"
        )

    # La identidad esta incluida.
    # Por ello utilizamos la proporcion exacta.
    p_exacto = np.mean(
        valores_F >= (F_obs - 1e-12)
    )

    print("\np exacto para tratamiento:")
    print(f"{p_exacto:.8f}")

    # ----------------------------------------
    # Guardar distribucion
    # ----------------------------------------

    archivo_perm = (
        "PERMANOVA_permutaciones_"
        + nombre_analisis
        + ".csv"
    )

    pd.DataFrame(
        {
            "F_permuted": valores_F
        }
    ).to_csv(
        archivo_perm,
        index=False
    )

    # ----------------------------------------
    # Resultado
    # ----------------------------------------

    fila = {
        "analysis": nombre_analisis,

        "n_signatures":
            Xdatos.shape[0],

        "n_features":
            Xdatos.shape[1],

        "distance":
            "Euclidean",

        "permutation_scheme":
            "exact within experimental block",

        "n_exact_permutations":
            contador,

        "F_treatment":
            observado["F"],

        "df_treatment":
            observado["df_effect"],

        "df_residual":
            observado["df_residual"],

        "SS_treatment":
            observado["SS_effect"],

        "SS_residual":
            observado["SS_residual"],

        "R2_treatment":
            observado["R2_total"],

        "partial_R2_treatment":
            observado["partial_R2"],

        "p_exact_treatment":
            p_exacto,

        "SS_block_context":
            contexto["SS_block_context"],

        "R2_block_context":
            contexto["R2_block_context"],

        "p_block_context":
            np.nan,
    }

    return fila


# ------------------------------------------------------------
# 7. Analisis primario: 12,328 genes
# ------------------------------------------------------------

resultado_full = ejecutar_permanova(
    X,
    "full_12328_genes"
)

# ------------------------------------------------------------
# 8. Analisis de sensibilidad: genes landmark
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
        "Se esperaban 978 genes landmark."
    )

resultado_landmark = ejecutar_permanova(
    X_landmark,
    "landmark_978_genes"
)

# ------------------------------------------------------------
# 9. Tabla final
# ------------------------------------------------------------

resultados = pd.DataFrame(
    [
        resultado_full,
        resultado_landmark,
    ]
)

resultados.to_csv(
    "PERMANOVA_resultados_dataset_primario.csv",
    index=False
)

# ------------------------------------------------------------
# 10. Informe metodologico
# ------------------------------------------------------------

with open(
    "PERMANOVA_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "PERMANOVA PAREADA EXACTA\n"
    )

    f.write(
        "=" * 60 + "\n\n"
    )

    f.write(
        "Dataset primario:\n"
    )

    f.write(
        "12 firmas, 4 bloques, "
        "3 tratamientos por bloque.\n"
    )

    f.write(
        "Tiempo: 24 h.\n"
    )

    f.write(
        "Dosis: 10 uM.\n\n"
    )

    f.write(
        "Distancia: Euclidea.\n"
    )

    f.write(
        "Modelo: tratamiento + bloque.\n"
    )

    f.write(
        "El efecto de tratamiento se evaluo "
        "controlando por bloque.\n"
    )

    f.write(
        "Las etiquetas de tratamiento fueron "
        "permutadas exclusivamente dentro "
        "de cada bloque experimental.\n"
    )

    f.write(
        "Se enumeraron exhaustivamente las "
        "1296 permutaciones posibles "
        "(3!^4).\n\n"
    )

    f.write(
        "El efecto bloque/contexto se cuantifico "
        "mediante R2 pero no se asigno un "
        "p inferencial, dado que cada linea "
        "celular corresponde a un unico "
        "bloque experimental en este "
        "subconjunto primario.\n\n"
    )

    f.write(
        "Se repitio el analisis utilizando "
        "exclusivamente los 978 genes "
        "landmark como analisis de robustez.\n"
    )

print("\n" + "=" * 72)

print("Archivos creados:")

print(
    " - PERMANOVA_resultados_dataset_primario.csv"
)

print(
    " - PERMANOVA_permutaciones_full_12328_genes.csv"
)

print(
    " - PERMANOVA_permutaciones_landmark_978_genes.csv"
)

print(
    " - PERMANOVA_metodologia.txt"
)

print("\nIMPORTANTE:")
print(
    "El p de tratamiento proviene de una "
    "permutacion EXACTA restringida por bloque."
)

print(
    "El efecto bloque/contexto se informa "
    "descriptivamente mediante R2."
)

print(
    "No se asigno p al efecto de linea celular "
    "en este subconjunto."
)

print(
    "\nPERMANOVA FINALIZADA CORRECTAMENTE"
)