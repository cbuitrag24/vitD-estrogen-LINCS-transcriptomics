import pandas as pd

# ============================================================
# CONSTRUCCION DEL DISEÑO PAREADO
# 24 h - 10 uM - LINCS L1000
#
# Este script usa SOLO metadatos.
# No inspecciona valores de expresion.
# ============================================================

ARCHIVO = "candidato_dataset_primario_10uM.csv"

d = pd.read_csv(ARCHIVO)

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol",
]

print("=" * 72)
print("CONSTRUCCION DEL DISEÑO PAREADO")
print("=" * 72)

# ------------------------------------------------------------
# 1. Obtener prefijo experimental
# ------------------------------------------------------------

def obtener_prefijo(sig_id):
    return str(sig_id).split("_")[0]


d["experiment_prefix"] = (
    d["sig_id"]
    .astype(str)
    .apply(obtener_prefijo)
)

# ------------------------------------------------------------
# 2. Definir bloque experimental
# ------------------------------------------------------------
# Un bloque queda definido por:
# linea celular + prefijo experimental

d["block_id"] = (
    d["cell_id"].astype(str)
    + "__"
    + d["experiment_prefix"].astype(str)
)

# ------------------------------------------------------------
# 3. Contar firmas por bloque y tratamiento
# ------------------------------------------------------------

conteo = (
    d.groupby(
        ["block_id", "cell_id",
         "experiment_prefix", "pert_iname"]
    )
    .size()
    .reset_index(name="n_firmas")
)

tabla = (
    conteo.pivot_table(
        index=[
            "block_id",
            "cell_id",
            "experiment_prefix"
        ],
        columns="pert_iname",
        values="n_firmas",
        fill_value=0
    )
    .reset_index()
)

# Garantizar que existan las tres columnas
for tratamiento in TRATAMIENTOS:
    if tratamiento not in tabla.columns:
        tabla[tratamiento] = 0

# ------------------------------------------------------------
# 4. Identificar bloques completos
# ------------------------------------------------------------

tabla["n_tratamientos"] = (
    (tabla["estradiol"] > 0).astype(int)
    + (tabla["calcitriol"] > 0).astype(int)
    + (tabla["calcipotriol"] > 0).astype(int)
)

tabla["triplete_completo"] = (
    (tabla["estradiol"] > 0)
    & (tabla["calcitriol"] > 0)
    & (tabla["calcipotriol"] > 0)
)

# Comprobar si existe exactamente una firma
# por tratamiento dentro del bloque

tabla["una_firma_por_tratamiento"] = (
    (tabla["estradiol"] == 1)
    & (tabla["calcitriol"] == 1)
    & (tabla["calcipotriol"] == 1)
)

# ------------------------------------------------------------
# 5. Comparaciones pareadas de dos tratamientos
# ------------------------------------------------------------

tabla["pareado_E2_calcitriol"] = (
    (tabla["estradiol"] > 0)
    & (tabla["calcitriol"] > 0)
)

tabla["pareado_E2_calcipotriol"] = (
    (tabla["estradiol"] > 0)
    & (tabla["calcipotriol"] > 0)
)

tabla["pareado_calcitriol_calcipotriol"] = (
    (tabla["calcitriol"] > 0)
    & (tabla["calcipotriol"] > 0)
)

# ------------------------------------------------------------
# 6. Guardar auditoria completa de bloques
# ------------------------------------------------------------

tabla.to_csv(
    "diseno_pareado_bloques.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Extraer bloques con los tres tratamientos
# ------------------------------------------------------------

bloques_triplete = tabla.loc[
    tabla["triplete_completo"],
    "block_id"
].tolist()

tripletes = d[
    d["block_id"].isin(bloques_triplete)
].copy()

tripletes = tripletes.sort_values(
    [
        "cell_id",
        "experiment_prefix",
        "pert_iname",
        "sig_id"
    ]
)

tripletes.to_csv(
    "candidato_primario_tripletes_10uM.csv",
    index=False
)

# ------------------------------------------------------------
# 8. Crear tabla compacta del candidato primario
# ------------------------------------------------------------

columnas = [
    "block_id",
    "cell_id",
    "experiment_prefix",
    "sig_id",
    "pert_id",
    "pert_iname",
    "pert_dose",
    "pert_dose_unit",
    "pert_time",
    "pert_time_unit",
    "distil_id",
]

columnas = [
    c for c in columnas
    if c in tripletes.columns
]

tripletes[columnas].to_csv(
    "candidato_primario_tripletes_detalle.csv",
    index=False
)

# ------------------------------------------------------------
# 9. Resumen de bloques completos
# ------------------------------------------------------------

resumen_tripletes = tabla[
    tabla["triplete_completo"]
].copy()

resumen_tripletes.to_csv(
    "candidato_primario_tripletes_resumen.csv",
    index=False
)

# ------------------------------------------------------------
# 10. Mostrar resultados
# ------------------------------------------------------------

print("\nNumero total de firmas candidatas a 10 uM:")
print(len(d))

print("\nNumero total de bloques experimentales:")
print(tabla["block_id"].nunique())

print("\nBloques con los TRES tratamientos:")
print(int(tabla["triplete_completo"].sum()))

print("\nDetalle de los bloques completos:")
print(
    tabla.loc[
        tabla["triplete_completo"],
        [
            "block_id",
            "cell_id",
            "experiment_prefix",
            "calcipotriol",
            "calcitriol",
            "estradiol",
            "una_firma_por_tratamiento"
        ]
    ].to_string(index=False)
)

print("\nFirmas dentro de los bloques completos:")
print(len(tripletes))

print("\nFirmas por tratamiento dentro de bloques completos:")
print(
    tripletes["pert_iname"]
    .value_counts()
    .to_string()
)

print("\nBloques pareados E2 vs calcitriol:")
print(int(tabla["pareado_E2_calcitriol"].sum()))

print("\nBloques pareados E2 vs calcipotriol:")
print(int(tabla["pareado_E2_calcipotriol"].sum()))

print("\nBloques pareados calcitriol vs calcipotriol:")
print(int(tabla["pareado_calcitriol_calcipotriol"].sum()))

print("\nArchivos creados:")
print(" - diseno_pareado_bloques.csv")
print(" - candidato_primario_tripletes_10uM.csv")
print(" - candidato_primario_tripletes_detalle.csv")
print(" - candidato_primario_tripletes_resumen.csv")

print("\nIMPORTANTE:")
print("La construccion de estos bloques se realizo")
print("exclusivamente con metadatos experimentales.")
print("No se inspeccionaron valores de expresion.")

print("\nDISEÑO PAREADO CONSTRUIDO CORRECTAMENTE")