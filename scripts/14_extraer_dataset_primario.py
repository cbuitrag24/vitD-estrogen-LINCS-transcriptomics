import numpy as np
import pandas as pd

# ============================================================
# EXTRACCION DEL DATASET PRIMARIO PAREADO
# LINCS L1000 - 24 h - 10 uM
#
# Este script:
# 1. Lee la matriz Level 5 previamente extraida.
# 2. Lee los metadatos correspondientes a las 213 firmas.
# 3. Lee las 12 firmas seleccionadas ANTES de inspeccionar
#    sus valores de expresion.
# 4. Extrae esas 12 firmas sin modificar sus valores.
# ============================================================

ARCHIVO_EXPRESION = "expresion_level5_213x12328.npz"
ARCHIVO_METADATA = "metadata_firmas_level5.csv"
ARCHIVO_GENES = "metadata_genes_level5.csv"
ARCHIVO_SELECCION = "candidato_primario_tripletes_10uM.csv"

print("=" * 70)
print("EXTRACCION DEL DATASET PRIMARIO PAREADO")
print("=" * 70)

# ------------------------------------------------------------
# 1. Cargar matriz de expresion
# ------------------------------------------------------------

npz = np.load(ARCHIVO_EXPRESION)

print("\nClaves disponibles en el archivo NPZ:")
print(npz.files)

# Detectar automaticamente la matriz
if "expression" in npz.files:
    X = npz["expression"]
elif "matrix" in npz.files:
    X = npz["matrix"]
elif "X" in npz.files:
    X = npz["X"]
else:
    # Si solo existe una matriz, usarla
    if len(npz.files) == 1:
        X = npz[npz.files[0]]
    else:
        raise ValueError(
            "No se pudo identificar automaticamente la matriz de expresion."
        )

print("\nForma de la matriz completa:")
print(X.shape)

# ------------------------------------------------------------
# 2. Cargar metadatos
# ------------------------------------------------------------

meta = pd.read_csv(ARCHIVO_METADATA)
genes = pd.read_csv(ARCHIVO_GENES)
seleccion = pd.read_csv(ARCHIVO_SELECCION)

print("\nFirmas en metadata:")
print(len(meta))

print("Genes en metadata:")
print(len(genes))

print("Firmas solicitadas para el dataset primario:")
print(len(seleccion))

# ------------------------------------------------------------
# 3. Controles de integridad
# ------------------------------------------------------------

if "sig_id" not in meta.columns:
    raise ValueError("metadata_firmas_level5.csv no contiene sig_id")

if "sig_id" not in seleccion.columns:
    raise ValueError("candidato_primario_tripletes_10uM.csv no contiene sig_id")

if meta["sig_id"].duplicated().any():
    raise ValueError("Hay sig_id duplicados en metadata_firmas_level5.csv")

if seleccion["sig_id"].duplicated().any():
    raise ValueError(
        "Hay sig_id duplicados en candidato_primario_tripletes_10uM.csv"
    )

faltantes = seleccion.loc[
    ~seleccion["sig_id"].isin(meta["sig_id"]),
    "sig_id"
]

print("\nFirmas seleccionadas ausentes de la matriz:")
print(len(faltantes))

if len(faltantes) > 0:
    print(faltantes.to_string(index=False))
    raise ValueError(
        "Hay firmas del diseño primario que no están en la matriz Level 5."
    )

# Verificar orientacion esperada
if X.shape[0] != len(meta):
    raise ValueError(
        "El numero de filas de la matriz no coincide con "
        "metadata_firmas_level5.csv"
    )

if X.shape[1] != len(genes):
    raise ValueError(
        "El numero de columnas de la matriz no coincide con "
        "metadata_genes_level5.csv"
    )

# ------------------------------------------------------------
# 4. Obtener posiciones de las 12 firmas
# ------------------------------------------------------------

mapa_posiciones = pd.Series(
    np.arange(len(meta)),
    index=meta["sig_id"]
)

indices = mapa_posiciones.loc[seleccion["sig_id"]].to_numpy()

# ------------------------------------------------------------
# 5. Extraer expresion SIN transformar los valores
# ------------------------------------------------------------

X_primario = X[indices, :].copy()

print("\nForma del dataset primario:")
print(X_primario.shape)

# ------------------------------------------------------------
# 6. Construir metadata en exactamente el mismo orden
# ------------------------------------------------------------

meta_primario = (
    seleccion[["sig_id"]]
    .merge(meta, on="sig_id", how="left", validate="one_to_one")
)

# Conservar columnas adicionales de la seleccion que no esten
# ya presentes en metadata

columnas_extra = [
    c for c in seleccion.columns
    if c != "sig_id" and c not in meta_primario.columns
]

if columnas_extra:
    extras = seleccion[["sig_id"] + columnas_extra]

    meta_primario = meta_primario.merge(
        extras,
        on="sig_id",
        how="left",
        validate="one_to_one"
    )

# ------------------------------------------------------------
# 7. Ultimos controles
# ------------------------------------------------------------

assert X_primario.shape[0] == len(meta_primario)
assert X_primario.shape[1] == len(genes)
assert list(meta_primario["sig_id"]) == list(seleccion["sig_id"])

# ------------------------------------------------------------
# 8. Guardar
# ------------------------------------------------------------

np.savez_compressed(
    "dataset_primario_expresion_12x12328.npz",
    expresion=X_primario
)

meta_primario.to_csv(
    "dataset_primario_metadata.csv",
    index=False
)

genes.to_csv(
    "dataset_primario_genes.csv",
    index=False
)

# ------------------------------------------------------------
# 9. Resumen
# ------------------------------------------------------------

print("\nFirmas por tratamiento:")

if "pert_iname" in meta_primario.columns:
    print(meta_primario["pert_iname"].value_counts())
else:
    print("No se encontro la columna pert_iname.")

print("\nLineas celulares:")

if "cell_id" in meta_primario.columns:
    print(sorted(meta_primario["cell_id"].unique()))

print("\nArchivos creados:")
print(" - dataset_primario_expresion_12x12328.npz")
print(" - dataset_primario_metadata.csv")
print(" - dataset_primario_genes.csv")

print("\nIMPORTANTE:")
print("Los valores de expresion fueron extraidos sin transformacion.")
print("No se seleccionaron firmas utilizando sus valores de expresion.")

print("\nEXTRACCION DEL DATASET PRIMARIO FINALIZADA CORRECTAMENTE")