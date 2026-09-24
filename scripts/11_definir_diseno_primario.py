import pandas as pd

# ============================================================
# DEFINICION A PRIORI DEL DISEÑO COMPARATIVO
# Estradiol vs vitamina D - LINCS L1000
# ============================================================

firmas = pd.read_csv("metadata_firmas_level5.csv")

# Siete lineas celulares que contienen:
# estradiol + calcitriol + calcipotriol
CELULAS_COMUNES = [
    "A375",
    "A549",
    "HA1E",
    "HT29",
    "MCF7",
    "PC3",
    "VCAP",
]

# ------------------------------------------------------------
# 1. Dataset de las siete lineas comunes
# ------------------------------------------------------------

comunes = firmas[
    firmas["cell_id"].isin(CELULAS_COMUNES)
].copy()

comunes["pert_dose"] = pd.to_numeric(
    comunes["pert_dose"],
    errors="coerce"
)

# Orden reproducible
comunes = comunes.sort_values(
    ["cell_id", "pert_iname", "pert_dose", "sig_id"]
)

comunes.to_csv(
    "diseno_7_celulas_todas_dosis.csv",
    index=False
)

# ------------------------------------------------------------
# 2. Tabla detallada de cada firma
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
    "distil_id",
]

columnas = [
    c for c in columnas
    if c in comunes.columns
]

detalle = comunes[columnas].copy()

detalle.to_csv(
    "diseno_detallado_7_celulas.csv",
    index=False
)

# ------------------------------------------------------------
# 3. Conteo por celula, compuesto y dosis
# ------------------------------------------------------------

conteo = (
    comunes
    .groupby(
        ["cell_id", "pert_iname", "pert_dose"],
        dropna=False
    )
    .size()
    .reset_index(name="n_firmas")
)

conteo.to_csv(
    "diseno_conteo_por_dosis.csv",
    index=False
)

# ------------------------------------------------------------
# 4. Identificar dosis compartidas dentro de cada celula
# ------------------------------------------------------------

registros = []

for celula in CELULAS_COMUNES:

    sub = comunes[
        comunes["cell_id"] == celula
    ]

    dosis_e2 = set(
        sub.loc[
            sub["pert_iname"] == "estradiol",
            "pert_dose"
        ].dropna()
    )

    dosis_calcitriol = set(
        sub.loc[
            sub["pert_iname"] == "calcitriol",
            "pert_dose"
        ].dropna()
    )

    dosis_calcipotriol = set(
        sub.loc[
            sub["pert_iname"] == "calcipotriol",
            "pert_dose"
        ].dropna()
    )

    compartidas_3 = (
        dosis_e2
        & dosis_calcitriol
        & dosis_calcipotriol
    )

    registros.append(
        {
            "cell_id": celula,
            "dosis_E2": "; ".join(
                map(str, sorted(dosis_e2))
            ),
            "dosis_calcitriol": "; ".join(
                map(str, sorted(dosis_calcitriol))
            ),
            "dosis_calcipotriol": "; ".join(
                map(str, sorted(dosis_calcipotriol))
            ),
            "dosis_compartidas_3": "; ".join(
                map(str, sorted(compartidas_3))
            ),
        }
    )

comparacion = pd.DataFrame(registros)

comparacion.to_csv(
    "diseno_dosis_compartidas.csv",
    index=False
)

# ------------------------------------------------------------
# 5. Resumen
# ------------------------------------------------------------

print("=" * 65)
print("DISEÑO PREVIO AL ANALISIS TRANSCRIPTOMICO")
print("=" * 65)

print("\nLineas celulares incluidas:")
print(", ".join(CELULAS_COMUNES))

print("\nNumero de firmas en las 7 lineas:")
print(len(comunes))

print("\nFirmas por compuesto:")
print(
    comunes["pert_iname"]
    .value_counts()
    .to_string()
)

print("\nFirmas por celula y compuesto:")
print(
    pd.crosstab(
        comunes["cell_id"],
        comunes["pert_iname"]
    ).to_string()
)

print("\nDosis compartidas por los tres tratamientos:")
print(
    comparacion[
        ["cell_id", "dosis_compartidas_3"]
    ].to_string(index=False)
)

print("\nArchivos creados:")
print(" - diseno_7_celulas_todas_dosis.csv")
print(" - diseno_detallado_7_celulas.csv")
print(" - diseno_conteo_por_dosis.csv")
print(" - diseno_dosis_compartidas.csv")

print("\nIMPORTANTE:")
print("Todavia no se utilizaron los valores de expresion")
print("para seleccionar firmas o dosis.")