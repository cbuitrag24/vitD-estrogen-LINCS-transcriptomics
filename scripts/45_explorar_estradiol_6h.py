import pandas as pd
from pathlib import Path

print("=" * 90)
print("EXPLORACIÓN DE ESTRADIOL A 6 H EN LINCS L1000")
print("=" * 90)

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

print("\nTotal de firmas en sig_info:", len(df))

# ------------------------------------------------------------
# 3. Seleccionar estradiol canónico
# ------------------------------------------------------------

# BRD-K18910433 es el pert_id de estradiol canónico
# que ya fue utilizado en nuestra curación previa.

e2 = df[
    (df["pert_id"] == "BRD-K18910433") &
    (df["pert_iname"].str.lower() == "estradiol")
].copy()

print("\nFirmas de estradiol canónico:", len(e2))

# ------------------------------------------------------------
# 4. Seleccionar 6 h
# ------------------------------------------------------------

e2_6h = e2[
    (pd.to_numeric(e2["pert_time"], errors="coerce") == 6) &
    (e2["pert_time_unit"].str.lower() == "h")
].copy()

print("\nFirmas de estradiol a 6 h:", len(e2_6h))

if len(e2_6h) == 0:
    print("\nNO HAY FIRMAS DE ESTRADIOL A 6 H.")
    print("No se justifica análisis temporal con este dataset.")
    raise SystemExit

# ------------------------------------------------------------
# 5. Número de perfiles constituyentes
# ------------------------------------------------------------

def contar_distil(x):
    if pd.isna(x):
        return 0

    x = str(x).strip()

    if not x:
        return 0

    return x.count("|") + 1


e2_6h["n_distil"] = e2_6h["distil_id"].apply(
    contar_distil
)

e2_6h["cumple_min3"] = e2_6h["n_distil"] >= 3

# ------------------------------------------------------------
# 6. Resumen por célula
# ------------------------------------------------------------

print("\n" + "=" * 90)
print("RESUMEN POR CÉLULA")
print("=" * 90)

resumen_celula = (
    e2_6h
    .groupby("cell_id")
    .agg(
        n_firmas=("sig_id", "count"),
        n_firmas_min3=("cumple_min3", "sum"),
        min_distil=("n_distil", "min"),
        max_distil=("n_distil", "max")
    )
    .reset_index()
    .sort_values(
        ["n_firmas_min3", "n_firmas"],
        ascending=False
    )
)

print(
    resumen_celula.to_string(index=False)
)

# ------------------------------------------------------------
# 7. Resumen por célula y dosis
# ------------------------------------------------------------

print("\n" + "=" * 90)
print("RESUMEN POR CÉLULA Y DOSIS")
print("=" * 90)

e2_6h["dose_numeric"] = pd.to_numeric(
    e2_6h["pert_dose"],
    errors="coerce"
)

resumen_dosis = (
    e2_6h
    .groupby(
        [
            "cell_id",
            "dose_numeric",
            "pert_dose_unit"
        ],
        dropna=False
    )
    .agg(
        n_firmas=("sig_id", "count"),
        n_firmas_min3=("cumple_min3", "sum")
    )
    .reset_index()
    .sort_values(
        ["cell_id", "dose_numeric"]
    )
)

print(
    resumen_dosis.to_string(index=False)
)

# ------------------------------------------------------------
# 8. Contextos del análisis primario
# ------------------------------------------------------------

primarias = ["A549", "HA1E", "PC3", "VCAP"]

print("\n" + "=" * 90)
print("COBERTURA EN LOS 4 CONTEXTOS DEL ANÁLISIS PRIMARIO")
print("=" * 90)

sub = e2_6h[
    e2_6h["cell_id"].isin(primarias)
].copy()

for cell in primarias:

    temp = sub[
        sub["cell_id"] == cell
    ]

    print(f"\n{cell}")
    print("-" * 40)

    if len(temp) == 0:
        print("Sin firmas E2 a 6 h.")
        continue

    columnas = [
        "sig_id",
        "pert_dose",
        "pert_dose_unit",
        "n_distil",
        "cumple_min3"
    ]

    print(
        temp[columnas].to_string(index=False)
    )

# ------------------------------------------------------------
# 9. Cobertura estricta >=3
# ------------------------------------------------------------

validas = e2_6h[
    e2_6h["cumple_min3"]
].copy()

print("\n" + "=" * 90)
print("RESUMEN FINAL")
print("=" * 90)

print(
    "Firmas E2 6 h con >=3 perfiles constituyentes:",
    len(validas)
)

print(
    "Número de células con al menos una firma válida:",
    validas["cell_id"].nunique()
)

print(
    "Células:",
    ", ".join(
        sorted(validas["cell_id"].unique())
    )
)

primarias_validas = sorted(
    set(validas["cell_id"]) &
    set(primarias)
)

print(
    "Contextos primarios con cobertura válida:",
    ", ".join(primarias_validas)
    if primarias_validas
    else "ninguno"
)

print(
    "Número de contextos primarios cubiertos:",
    len(primarias_validas),
    "/ 4"
)

# ------------------------------------------------------------
# 10. Guardar resultados
# ------------------------------------------------------------

e2_6h.to_csv(
    "ESTRADIOL_6H_firmas.csv",
    index=False
)

resumen_celula.to_csv(
    "ESTRADIOL_6H_resumen_celulas.csv",
    index=False
)

resumen_dosis.to_csv(
    "ESTRADIOL_6H_resumen_dosis.csv",
    index=False
)

print("\nArchivos generados:")
print("- ESTRADIOL_6H_firmas.csv")
print("- ESTRADIOL_6H_resumen_celulas.csv")
print("- ESTRADIOL_6H_resumen_dosis.csv")

print("\nFIN")