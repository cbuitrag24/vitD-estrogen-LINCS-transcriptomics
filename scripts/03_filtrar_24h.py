import pandas as pd

# --------------------------------------------------
# 1. Cargar SIG_INFO
# --------------------------------------------------

sig_file = "data/GSE92742_Broad_LINCS_sig_info.txt.gz"

sig = pd.read_csv(
    sig_file,
    sep="\t",
    compression="gzip",
    low_memory=False
)

# --------------------------------------------------
# 2. Definir los compuestos principales
# --------------------------------------------------

compuestos = [
    "estradiol",
    "17-beta-estradiol",
    "calcitriol",
    "calcipotriol"
]

# --------------------------------------------------
# 3. Quedarnos solo con esos compuestos
# --------------------------------------------------

datos = sig[
    sig["pert_iname"]
    .astype(str)
    .str.lower()
    .isin(compuestos)
].copy()

# --------------------------------------------------
# 4. Filtrar solo tratamientos de 24 horas
# --------------------------------------------------

datos_24h = datos[
    (datos["pert_time"] == 24) &
    (datos["pert_time_unit"].astype(str).str.lower() == "h")
].copy()

# --------------------------------------------------
# 5. Contar replicas desde distil_id
# --------------------------------------------------

def contar_replicas(valor):
    if pd.isna(valor):
        return 0

    texto = str(valor)

    # Los identificadores de replicas suelen estar separados
    # por el caracter |
    if "|" in texto:
        return len(texto.split("|"))

    # En algunos casos puede haber otro separador
    if "," in texto:
        return len(texto.split(","))

    return 1


datos_24h["n_replicas"] = datos_24h["distil_id"].apply(contar_replicas)

# --------------------------------------------------
# 6. Mostrar cantidad total de firmas a 24 h
# --------------------------------------------------

print("\n===================================")
print("FIRMAS A 24 HORAS")
print("===================================")

print(
    datos_24h
    .groupby("pert_iname")
    .size()
    .sort_values(ascending=False)
)

# --------------------------------------------------
# 7. Mostrar lineas celulares por compuesto
# --------------------------------------------------

print("\n===================================")
print("LINEAS CELULARES A 24 HORAS")
print("===================================")

lineas = (
    datos_24h[
        ["pert_iname", "cell_id"]
    ]
    .drop_duplicates()
    .sort_values(["pert_iname", "cell_id"])
)

print(lineas.to_string(index=False))

# --------------------------------------------------
# 8. Mostrar distribucion de replicas
# --------------------------------------------------

print("\n===================================")
print("NUMERO DE REPLICAS POR FIRMA")
print("===================================")

print(
    datos_24h[
        [
            "pert_iname",
            "pert_id",
            "cell_id",
            "pert_dose",
            "pert_dose_unit",
            "n_replicas"
        ]
    ]
    .sort_values(
        ["pert_iname", "cell_id", "pert_dose"]
    )
    .to_string(index=False)
)

# --------------------------------------------------
# 9. Filtrar firmas con al menos 3 replicas
# --------------------------------------------------

datos_24h_rep3 = datos_24h[
    datos_24h["n_replicas"] >= 3
].copy()

print("\n===================================")
print("FIRMAS 24 H CON >= 3 REPLICAS")
print("===================================")

print(
    datos_24h_rep3
    .groupby("pert_iname")
    .size()
    .sort_values(ascending=False)
)

# --------------------------------------------------
# 10. Lineas celulares despues del filtro
# --------------------------------------------------

print("\n===================================")
print("LINEAS CELULARES CON >= 3 REPLICAS")
print("===================================")

lineas_rep3 = (
    datos_24h_rep3[
        ["pert_iname", "cell_id"]
    ]
    .drop_duplicates()
    .sort_values(["pert_iname", "cell_id"])
)

print(lineas_rep3.to_string(index=False))

# --------------------------------------------------
# 11. Guardar tabla filtrada
# --------------------------------------------------

datos_24h_rep3.to_csv(
    "firmas_24h_rep3.csv",
    index=False
)

print("\nArchivo creado: firmas_24h_rep3.csv")