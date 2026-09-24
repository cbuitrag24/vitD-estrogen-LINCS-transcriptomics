import pandas as pd
import re

# ============================================================
# AUDITORIA DE REPLICACION
# Dataset candidato primario: 10 uM, 24 h
# ============================================================

ARCHIVO = "candidato_dataset_primario_10uM.csv"

d = pd.read_csv(ARCHIVO)

print("=" * 70)
print("AUDITORIA DE REPLICACION - DATASET CANDIDATO PRIMARIO")
print("=" * 70)

print("\nNumero total de firmas:", len(d))
print("Numero de sig_id unicos:", d["sig_id"].nunique())

# ------------------------------------------------------------
# 1. Contar perfiles constituyentes de cada firma Level 5
# ------------------------------------------------------------

def contar_distil(x):
    if pd.isna(x):
        return 0

    texto = str(x).strip()

    if texto == "":
        return 0

    # En LINCS los distil_id pueden estar separados por |
    return len(texto.split("|"))


d["n_distil"] = d["distil_id"].apply(contar_distil)

# ------------------------------------------------------------
# 2. Extraer prefijo experimental del sig_id
# ------------------------------------------------------------

def extraer_prefijo(sig):
    sig = str(sig)

    # Parte anterior al primer "_"
    # Ejemplo:
    # CPC011_A549_24H:... -> CPC011
    # MUC.CP001_MCF7_24H:... -> MUC.CP001

    return sig.split("_")[0]


d["experiment_prefix"] = d["sig_id"].apply(extraer_prefijo)

# ------------------------------------------------------------
# 3. Tabla firma por firma
# ------------------------------------------------------------

columnas = [
    "sig_id",
    "pert_id",
    "pert_iname",
    "cell_id",
    "pert_dose",
    "pert_dose_unit",
    "pert_time",
    "pert_time_unit",
    "experiment_prefix",
    "n_distil",
    "distil_id",
]

columnas = [c for c in columnas if c in d.columns]

d[columnas].to_csv(
    "auditoria_firma_por_firma.csv",
    index=False
)

# ------------------------------------------------------------
# 4. Resumen por celula y tratamiento
# ------------------------------------------------------------

resumen = (
    d.groupby(["cell_id", "pert_iname"])
    .agg(
        n_firmas=("sig_id", "nunique"),
        n_experimentos=("experiment_prefix", "nunique"),
        min_distil=("n_distil", "min"),
        mediana_distil=("n_distil", "median"),
        max_distil=("n_distil", "max"),
    )
    .reset_index()
)

resumen.to_csv(
    "auditoria_resumen_replicacion.csv",
    index=False
)

# ------------------------------------------------------------
# 5. Prefijos experimentales por celula y tratamiento
# ------------------------------------------------------------

prefijos = (
    d.groupby(["cell_id", "pert_iname"])["experiment_prefix"]
    .apply(
        lambda x: "; ".join(
            sorted(set(map(str, x)))
        )
    )
    .reset_index(name="experiment_prefixes")
)

prefijos.to_csv(
    "auditoria_prefijos_experimentales.csv",
    index=False
)

# ------------------------------------------------------------
# 6. Comprobar si un mismo prefijo contiene varios tratamientos
# ------------------------------------------------------------

cruce = (
    d.groupby(["cell_id", "experiment_prefix"])["pert_iname"]
    .agg(
        lambda x: "; ".join(
            sorted(set(map(str, x)))
        )
    )
    .reset_index(name="tratamientos")
)

cruce["n_tratamientos"] = (
    cruce["tratamientos"]
    .str.split("; ")
    .str.len()
)

cruce.to_csv(
    "auditoria_experimentos_tratamientos.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Mostrar resultados
# ------------------------------------------------------------

print("\nResumen por linea celular y tratamiento:")
print(resumen.to_string(index=False))

print("\nNumero de firmas segun cantidad de perfiles distil:")
print(
    d["n_distil"]
    .value_counts()
    .sort_index()
    .to_string()
)

print("\nExperimentos que contienen mas de un tratamiento:")

multi = cruce[
    cruce["n_tratamientos"] > 1
]

if len(multi) == 0:
    print("NINGUNO")
else:
    print(multi.to_string(index=False))

print("\nArchivos creados:")
print(" - auditoria_firma_por_firma.csv")
print(" - auditoria_resumen_replicacion.csv")
print(" - auditoria_prefijos_experimentales.csv")
print(" - auditoria_experimentos_tratamientos.csv")

print("\nAUDITORIA FINALIZADA CORRECTAMENTE")