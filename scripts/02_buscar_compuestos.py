import pandas as pd

# --------------------------------------------------
# 1. Cargar los archivos de metadata
# --------------------------------------------------

sig_file = "data/GSE92742_Broad_LINCS_sig_info.txt.gz"
pert_file = "data/GSE92742_Broad_LINCS_pert_info.txt.gz"

sig = pd.read_csv(
    sig_file,
    sep="\t",
    compression="gzip",
    low_memory=False
)

pert = pd.read_csv(
    pert_file,
    sep="\t",
    compression="gzip",
    low_memory=False
)

# --------------------------------------------------
# 2. Definir los compuestos que queremos buscar
# --------------------------------------------------

estrogenos = [
    "estradiol",
    "17-beta-estradiol",
    "alpha-estradiol",
    "ethinyl-estradiol",
    "estradiol-benzoate",
    "estradiol-valerate",
    "estradiol-cypionate"
]

vitamina_d = [
    "calcitriol",
    "calcipotriol",
    "maxacalcitol",
    "seocalcitol",
    "ercalcitriol",
    "tacalcitol",
    "paricalcitol"
]

# --------------------------------------------------
# 3. Buscar los compuestos en PERT_INFO
# --------------------------------------------------

estrogenos_pert = pert[
    pert["pert_iname"].astype(str).str.lower().isin(estrogenos)
].copy()

vitd_pert = pert[
    pert["pert_iname"].astype(str).str.lower().isin(vitamina_d)
].copy()

print("\n==============================")
print("ESTROGENOS EN PERT_INFO")
print("==============================")

print(
    estrogenos_pert[
        ["pert_id", "pert_iname", "pert_type"]
    ].to_string(index=False)
)

print("\n==============================")
print("VITAMINA D EN PERT_INFO")
print("==============================")

print(
    vitd_pert[
        ["pert_id", "pert_iname", "pert_type"]
    ].to_string(index=False)
)

# --------------------------------------------------
# 4. Buscar las firmas correspondientes en SIG_INFO
# --------------------------------------------------

todos_compuestos = estrogenos + vitamina_d

firmas = sig[
    sig["pert_iname"].astype(str).str.lower().isin(todos_compuestos)
].copy()

print("\n==============================")
print("NUMERO DE FIRMAS POR COMPUESTO")
print("==============================")

conteo = (
    firmas
    .groupby("pert_iname")
    .size()
    .sort_values(ascending=False)
)

print(conteo)

# --------------------------------------------------
# 5. Mostrar tiempos experimentales
# --------------------------------------------------

print("\n==============================")
print("TIEMPOS DISPONIBLES")
print("==============================")

tiempos = (
    firmas[
        ["pert_iname", "pert_time", "pert_time_unit"]
    ]
    .drop_duplicates()
    .sort_values(["pert_iname", "pert_time"])
)

print(tiempos.to_string(index=False))

# --------------------------------------------------
# 6. Mostrar lineas celulares
# --------------------------------------------------

print("\n==============================")
print("LINEAS CELULARES")
print("==============================")

celulas = (
    firmas[
        ["pert_iname", "cell_id"]
    ]
    .drop_duplicates()
    .sort_values(["pert_iname", "cell_id"])
)

print(celulas.to_string(index=False))