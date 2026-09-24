import pandas as pd

# Leer las firmas de 24 h con al menos 3 perfiles constituyentes
df = pd.read_csv("firmas_24h_rep3.csv")

# Estradiol canónico (E2)
e2 = df[df["pert_id"] == "BRD-K18910433"].copy()

# Firmas de vitamina D
vitd = df[df["pert_iname"].isin(["calcitriol", "calcipotriol"])].copy()

# Unir los grupos
seleccion = pd.concat([e2, vitd], ignore_index=True)

# Identificar la clase hormonal
seleccion["hormone_class"] = seleccion["pert_iname"].apply(
    lambda x: "E2" if x == "estradiol" else "Vitamin_D"
)

# Guardar el archivo
seleccion.to_csv("firmas_seleccionadas_level5.csv", index=False)

print("========================================")
print("SELECCION PARA LEVEL 5")
print("========================================")

print("\nFirmas por compuesto:")
print(seleccion["pert_iname"].value_counts())

print("\nFirmas por pert_id:")
print(
    seleccion.groupby(["pert_iname", "pert_id"])
    .size()
    .to_string()
)

print("\nTotal de firmas:", len(seleccion))
print("\nArchivo creado: firmas_seleccionadas_level5.csv")