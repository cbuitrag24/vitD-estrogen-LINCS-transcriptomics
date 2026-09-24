import pandas as pd
from pathlib import Path

print("=" * 95)
print("AUDITORÍA DEL DISEÑO TEMPORAL DE ESTRADIOL: 6 H vs 24 H")
print("=" * 95)

# ------------------------------------------------------------
# 1. Buscar sig_info
# ------------------------------------------------------------

posibles = [
    Path("GSE92742_Broad_LINCS_sig_info.txt.gz"),
    Path("data/GSE92742_Broad_LINCS_sig_info.txt.gz"),
]

archivo = None

for p in posibles:
    if p.exists():
        archivo = p
        break

if archivo is None:
    raise FileNotFoundError(
        "No se encontró GSE92742_Broad_LINCS_sig_info.txt.gz"
    )

print("\nArchivo utilizado:")
print(archivo)

# ------------------------------------------------------------
# 2. Leer metadatos
# ------------------------------------------------------------

df = pd.read_csv(
    archivo,
    sep="\t",
    low_memory=False
)

# ------------------------------------------------------------
# 3. Seleccionar estradiol canónico
# ------------------------------------------------------------

e2 = df[
    (df["pert_id"] == "BRD-K18910433") &
    (df["pert_iname"].str.lower() == "estradiol")
].copy()

e2["tiempo"] = pd.to_numeric(
    e2["pert_time"],
    errors="coerce"
)

e2["dosis"] = pd.to_numeric(
    e2["pert_dose"],
    errors="coerce"
)

# ------------------------------------------------------------
# 4. Contar perfiles constituyentes
# ------------------------------------------------------------

def contar_distil(x):
    if pd.isna(x):
        return 0

    x = str(x).strip()

    if not x:
        return 0

    return x.count("|") + 1


e2["n_distil"] = e2["distil_id"].apply(contar_distil)

# ------------------------------------------------------------
# 5. Selección estricta para sensibilidad temporal
# ------------------------------------------------------------

celulas = ["A549", "HA1E", "PC3", "VCAP"]

temporal = e2[
    (e2["cell_id"].isin(celulas)) &
    (e2["tiempo"].isin([6, 24])) &
    (e2["dosis"] == 10) &
    (e2["pert_dose_unit"] == "µM") &
    (e2["n_distil"] >= 3)
].copy()

# ------------------------------------------------------------
# 6. Obtener bloque experimental desde sig_id
# ------------------------------------------------------------

temporal["bloque"] = (
    temporal["sig_id"]
    .str.split("_")
    .str[0]
)

temporal = temporal.sort_values(
    ["cell_id", "bloque", "tiempo", "sig_id"]
)

print("\n" + "=" * 95)
print("FIRMAS VÁLIDAS: E2 10 µM, >=3 PERFILES")
print("=" * 95)

print(
    temporal[
        [
            "cell_id",
            "bloque",
            "tiempo",
            "sig_id",
            "n_distil"
        ]
    ].to_string(index=False)
)

# ------------------------------------------------------------
# 7. Auditar bloques con ambos tiempos
# ------------------------------------------------------------

conteo = (
    temporal
    .groupby(
        ["cell_id", "bloque", "tiempo"]
    )
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

if 6 not in conteo.columns:
    conteo[6] = 0

if 24 not in conteo.columns:
    conteo[24] = 0

conteo["tiene_6h"] = conteo[6] > 0
conteo["tiene_24h"] = conteo[24] > 0

conteo["pareado_6h_24h"] = (
    conteo["tiene_6h"] &
    conteo["tiene_24h"]
)

print("\n" + "=" * 95)
print("AUDITORÍA DE BLOQUES")
print("=" * 95)

print(
    conteo[
        [
            "cell_id",
            "bloque",
            6,
            24,
            "pareado_6h_24h"
        ]
    ].to_string(index=False)
)

# ------------------------------------------------------------
# 8. Bloques temporalmente apareados
# ------------------------------------------------------------

pareados = conteo[
    conteo["pareado_6h_24h"]
].copy()

print("\n" + "=" * 95)
print("BLOQUES CON E2 A 6 H Y 24 H")
print("=" * 95)

if len(pareados) == 0:

    print("No existen bloques temporalmente apareados.")

else:

    print(
        pareados[
            [
                "cell_id",
                "bloque",
                6,
                24
            ]
        ].to_string(index=False)
    )

# ------------------------------------------------------------
# 9. Firmas exactas de los bloques apareados
# ------------------------------------------------------------

if len(pareados) > 0:

    claves = pareados[
        ["cell_id", "bloque"]
    ].copy()

    seleccion = temporal.merge(
        claves,
        on=["cell_id", "bloque"],
        how="inner"
    )

    seleccion = seleccion.sort_values(
        ["cell_id", "bloque", "tiempo"]
    )

    print("\n" + "=" * 95)
    print("FIRMAS QUE PODRÍAN UTILIZARSE EN EL ANÁLISIS TEMPORAL")
    print("=" * 95)

    print(
        seleccion[
            [
                "cell_id",
                "bloque",
                "tiempo",
                "sig_id",
                "n_distil"
            ]
        ].to_string(index=False)
    )

    seleccion.to_csv(
        "E2_temporal_firmas_pareadas.csv",
        index=False
    )

# ------------------------------------------------------------
# 10. Resumen final
# ------------------------------------------------------------

print("\n" + "=" * 95)
print("RESUMEN FINAL")
print("=" * 95)

print(
    "Número de bloques temporalmente apareados:",
    len(pareados)
)

print(
    "Número de células representadas:",
    pareados["cell_id"].nunique()
)

if len(pareados) > 0:

    print(
        "Células:",
        ", ".join(
            sorted(pareados["cell_id"].unique())
        )
    )

    print(
        "Bloques:",
        ", ".join(
            sorted(pareados["bloque"].unique())
        )
    )

conteo.to_csv(
    "E2_temporal_auditoria_bloques.csv",
    index=False
)

print("\nArchivo generado:")
print("- E2_temporal_auditoria_bloques.csv")

if len(pareados) > 0:
    print("- E2_temporal_firmas_pareadas.csv")

print("\nFIN")