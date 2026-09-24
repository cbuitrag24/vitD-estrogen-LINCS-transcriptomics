import pandas as pd

# Leer el archivo que ya generamos
df = pd.read_csv("firmas_24h_rep3.csv")

# Nos quedamos solamente con los tres tratamientos principales
compuestos = ["estradiol", "calcitriol", "calcipotriol"]

df = df[df["pert_iname"].isin(compuestos)].copy()

# Líneas celulares presentes en los tres tratamientos
lineas_por_compuesto = {
    compuesto: set(df.loc[df["pert_iname"] == compuesto, "cell_id"])
    for compuesto in compuestos
}

lineas_comunes = set.intersection(*lineas_por_compuesto.values())

print("\n======================================")
print("LINEAS PRESENTES EN LOS TRES TRATAMIENTOS")
print("======================================")
print(sorted(lineas_comunes))

# Filtrar solamente esas líneas
comunes = df[df["cell_id"].isin(lineas_comunes)].copy()

print("\n======================================")
print("NUMERO DE FIRMAS")
print("======================================")

tabla_firmas = (
    comunes
    .groupby(["cell_id", "pert_iname"])
    .size()
    .unstack(fill_value=0)
)

print(tabla_firmas.to_string())

print("\n======================================")
print("DOSIS DISPONIBLES")
print("======================================")

for celula in sorted(lineas_comunes):
    print(f"\n--- {celula} ---")

    sub = comunes[comunes["cell_id"] == celula]

    for compuesto in compuestos:
        datos = sub[sub["pert_iname"] == compuesto]

        dosis = (
            datos[["pert_dose", "pert_dose_unit"]]
            .drop_duplicates()
            .sort_values("pert_dose")
        )

        print(f"{compuesto}:")
        if len(dosis) == 0:
            print("  sin datos")
        else:
            for _, fila in dosis.iterrows():
                print(
                    f"  {fila['pert_dose']} {fila['pert_dose_unit']}"
                )

# Guardar la tabla
tabla_firmas.to_csv("tabla_diseno_comparativo.csv")

print("\n======================================")
print("ARCHIVO CREADO")
print("======================================")
print("tabla_diseno_comparativo.csv")