import pandas as pd

# --------------------------------------------------
# 1. Cargar el archivo ya filtrado
# --------------------------------------------------

archivo = "firmas_24h_rep3.csv"

df = pd.read_csv(archivo)

print("\n===================================")
print("RESUMEN GENERAL")
print("===================================")

print("Numero total de firmas:", len(df))

print("\nFirmas por compuesto:")
print(
    df.groupby("pert_iname")
    .size()
    .sort_values(ascending=False)
)

# --------------------------------------------------
# 2. Ver los pert_id de cada compuesto
# --------------------------------------------------

print("\n===================================")
print("PERT_ID POR COMPUESTO")
print("===================================")

pert_ids = (
    df[
        ["pert_iname", "pert_id"]
    ]
    .drop_duplicates()
    .sort_values(["pert_iname", "pert_id"])
)

print(pert_ids.to_string(index=False))

# --------------------------------------------------
# 3. Contar firmas por compuesto y linea celular
# --------------------------------------------------

print("\n===================================")
print("FIRMAS POR COMPUESTO Y LINEA CELULAR")
print("===================================")

conteo_celulas = (
    df.groupby(
        ["pert_iname", "cell_id"]
    )
    .size()
    .reset_index(name="n_firmas")
    .sort_values(
        ["pert_iname", "n_firmas"],
        ascending=[True, False]
    )
)

print(conteo_celulas.to_string(index=False))

# --------------------------------------------------
# 4. Identificar lineas celulares por compuesto
# --------------------------------------------------

lineas_por_compuesto = {}

for compuesto in sorted(df["pert_iname"].unique()):

    lineas = set(
        df.loc[
            df["pert_iname"] == compuesto,
            "cell_id"
        ]
    )

    lineas_por_compuesto[compuesto] = lineas

# --------------------------------------------------
# 5. Mostrar intersecciones
# --------------------------------------------------

print("\n===================================")
print("LINEAS COMPARTIDAS")
print("===================================")

estradiol = lineas_por_compuesto.get("estradiol", set())
beta_e2 = lineas_por_compuesto.get("17-beta-estradiol", set())
calcitriol = lineas_por_compuesto.get("calcitriol", set())
calcipotriol = lineas_por_compuesto.get("calcipotriol", set())

print("\nEstradiol vs calcitriol:")
print(sorted(estradiol & calcitriol))

print("\nEstradiol vs calcipotriol:")
print(sorted(estradiol & calcipotriol))

print("\nEstradiol vs ambos vitamina D:")
print(
    sorted(
        estradiol &
        calcitriol &
        calcipotriol
    )
)

print("\n17-beta-estradiol vs calcitriol:")
print(sorted(beta_e2 & calcitriol))

print("\n17-beta-estradiol vs calcipotriol:")
print(sorted(beta_e2 & calcipotriol))

# --------------------------------------------------
# 6. Dosis disponibles
# --------------------------------------------------

print("\n===================================")
print("DOSIS DISPONIBLES")
print("===================================")

dosis = (
    df[
        [
            "pert_iname",
            "cell_id",
            "pert_dose",
            "pert_dose_unit"
        ]
    ]
    .drop_duplicates()
    .sort_values(
        ["pert_iname", "cell_id", "pert_dose"]
    )
)

print(dosis.to_string(index=False))

# --------------------------------------------------
# 7. Guardar resumen en CSV
# --------------------------------------------------

conteo_celulas.to_csv(
    "resumen_firmas_por_celula.csv",
    index=False
)

dosis.to_csv(
    "resumen_dosis.csv",
    index=False
)

print("\nArchivos creados:")
print("resumen_firmas_por_celula.csv")
print("resumen_dosis.csv")