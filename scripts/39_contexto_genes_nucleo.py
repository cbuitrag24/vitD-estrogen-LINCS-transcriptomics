import numpy as np
import pandas as pd

print("=" * 90)
print("ANÁLISIS POR CONTEXTO DE GENES NÚCLEO")
print("=" * 90)

# ================================================================
# 1. ARCHIVOS
# ================================================================

ARCHIVO_EXPR = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_GENES = "dataset_primario_genes.csv"
ARCHIVO_META = "dataset_primario_metadata.csv"

GENES_NUCLEO = [
    "CDC20",
    "MYC",
    "BIRC5",
    "CCNA2",
    "CCNB2",
    "CDC25A",
    "EZH2",
    "KIF2C",
    "LBR",
    "MCM3",
    "NOLC1",
    "PCNA",
    "USP1"
]

# ================================================================
# 2. CARGAR DATOS
# ================================================================

npz = np.load(ARCHIVO_EXPR)
X = npz["expresion"]

genes = pd.read_csv(ARCHIVO_GENES)
meta = pd.read_csv(ARCHIVO_META)

print("\nDimensión matriz de expresión:", X.shape)
print("Dimensión metadata:", meta.shape)
print("Dimensión tabla de genes:", genes.shape)

if X.shape[0] != len(meta):
    raise ValueError(
        f"No coincide número de firmas: matriz={X.shape[0]} metadata={len(meta)}"
    )

if X.shape[1] != len(genes):
    raise ValueError(
        f"No coincide número de genes: matriz={X.shape[1]} genes={len(genes)}"
    )

# ================================================================
# 3. IDENTIFICAR COLUMNAS
# ================================================================

if "pr_gene_symbol" not in genes.columns:
    raise ValueError("No encontré la columna pr_gene_symbol en dataset_primario_genes.csv")

if "cell_id" not in meta.columns:
    raise ValueError("No encontré cell_id en metadata")

if "tratamiento" in meta.columns:
    col_trat = "tratamiento"
elif "pert_iname" in meta.columns:
    col_trat = "pert_iname"
else:
    raise ValueError(
        "No encontré columna tratamiento ni pert_iname en metadata"
    )

print("\nColumna de tratamiento usada:", col_trat)
print("Tratamientos encontrados:", sorted(meta[col_trat].astype(str).unique()))
print("Contextos encontrados:", sorted(meta["cell_id"].astype(str).unique()))

# ================================================================
# 4. LOCALIZAR LOS GENES NÚCLEO
# ================================================================

simbolos = genes["pr_gene_symbol"].astype(str)

indices = {}
faltantes = []

for gene in GENES_NUCLEO:
    encontrados = np.where(simbolos.values == gene)[0]

    if len(encontrados) == 0:
        faltantes.append(gene)
    else:
        indices[gene] = encontrados[0]

print("\nGenes núcleo encontrados:", len(indices))
print(list(indices.keys()))

if faltantes:
    print("\nGenes no encontrados:")
    print(faltantes)

if len(indices) == 0:
    raise ValueError("No se encontró ninguno de los genes núcleo")

# ================================================================
# 5. EXTRAER VALORES INDIVIDUALES
# ================================================================

filas = []

for i in range(len(meta)):

    contexto = meta.iloc[i]["cell_id"]
    tratamiento = meta.iloc[i][col_trat]

    sig_id = (
        meta.iloc[i]["sig_id"]
        if "sig_id" in meta.columns
        else f"firma_{i}"
    )

    for gene, j in indices.items():

        filas.append({
            "sig_id": sig_id,
            "cell_id": contexto,
            "tratamiento": tratamiento,
            "gene": gene,
            "modz": float(X[i, j])
        })

valores = pd.DataFrame(filas)

valores.to_csv(
    "CONTEXTO_genes_nucleo_valores.csv",
    index=False
)

# ================================================================
# 6. NORMALIZAR NOMBRES DE TRATAMIENTO
# ================================================================

def normalizar_tratamiento(x):

    x = str(x).strip().lower()

    if "estradiol" in x:
        return "estradiol"

    if "calcipotriol" in x:
        return "calcipotriol"

    if "calcitriol" in x:
        return "calcitriol"

    return x


valores["tratamiento_normalizado"] = (
    valores["tratamiento"]
    .apply(normalizar_tratamiento)
)

print(
    "\nTratamientos normalizados:",
    sorted(valores["tratamiento_normalizado"].unique())
)

# ================================================================
# 7. RESUMEN POR CONTEXTO
# ================================================================

resumen = (
    valores
    .groupby(
        ["cell_id", "tratamiento_normalizado", "gene"],
        as_index=False
    )
    .agg(
        modz=("modz", "mean")
    )
)

resumen["signo"] = np.where(
    resumen["modz"] > 0,
    1,
    np.where(
        resumen["modz"] < 0,
        -1,
        0
    )
)

resumen.to_csv(
    "CONTEXTO_genes_nucleo_resumen.csv",
    index=False
)

# ================================================================
# 8. MATRIZ DE SIGNOS
# ================================================================

signos = (
    resumen
    .pivot_table(
        index=["gene", "tratamiento_normalizado"],
        columns="cell_id",
        values="signo",
        aggfunc="first"
    )
    .reset_index()
)

signos.to_csv(
    "CONTEXTO_genes_nucleo_signos.csv",
    index=False
)

# ================================================================
# 9. CONSISTENCIA ENTRE LOS 4 CONTEXTOS
# ================================================================

consistencia = []

for tratamiento in sorted(
    resumen["tratamiento_normalizado"].unique()
):

    sub_t = resumen[
        resumen["tratamiento_normalizado"] == tratamiento
    ]

    for gene in indices:

        sub_g = sub_t[
            sub_t["gene"] == gene
        ]

        valores_gene = sub_g["modz"].to_numpy()

        n = len(valores_gene)

        if n == 0:
            continue

        n_neg = int(np.sum(valores_gene < 0))
        n_pos = int(np.sum(valores_gene > 0))
        n_cero = int(np.sum(valores_gene == 0))

        if n_neg == n:
            clasificacion = f"negativo_{n}de{n}"

        elif n_pos == n:
            clasificacion = f"positivo_{n}de{n}"

        elif n_neg >= 3 and n == 4:
            clasificacion = "predominio_negativo_3de4"

        elif n_pos >= 3 and n == 4:
            clasificacion = "predominio_positivo_3de4"

        else:
            clasificacion = "heterogeneo"

        consistencia.append({
            "tratamiento": tratamiento,
            "gene": gene,
            "n_contextos": n,
            "n_negativos": n_neg,
            "n_positivos": n_pos,
            "n_cero": n_cero,
            "media_modz": float(np.mean(valores_gene)),
            "mediana_modz": float(np.median(valores_gene)),
            "min_modz": float(np.min(valores_gene)),
            "max_modz": float(np.max(valores_gene)),
            "clasificacion": clasificacion
        })

consistencia_df = pd.DataFrame(consistencia)

consistencia_df.to_csv(
    "CONTEXTO_genes_nucleo_consistencia.csv",
    index=False
)

# ================================================================
# 10. CONVERGENCIA ENTRE LOS 3 TRATAMIENTOS
# ================================================================

convergencia = []

tratamientos_objetivo = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

for gene in indices:

    sub = consistencia_df[
        consistencia_df["gene"] == gene
    ]

    presentes = set(sub["tratamiento"])

    if not set(tratamientos_objetivo).issubset(presentes):
        continue

    medias = dict(
        zip(
            sub["tratamiento"],
            sub["media_modz"]
        )
    )

    clasificaciones = dict(
        zip(
            sub["tratamiento"],
            sub["clasificacion"]
        )
    )

    vals = [
        medias[t]
        for t in tratamientos_objetivo
    ]

    if all(v < 0 for v in vals):
        direccion = "negativa_en_los_3"

    elif all(v > 0 for v in vals):
        direccion = "positiva_en_los_3"

    else:
        direccion = "mixta"

    convergencia.append({
        "gene": gene,
        "direccion_media": direccion,

        "estradiol_media":
            medias["estradiol"],

        "calcitriol_media":
            medias["calcitriol"],

        "calcipotriol_media":
            medias["calcipotriol"],

        "estradiol_clasificacion":
            clasificaciones["estradiol"],

        "calcitriol_clasificacion":
            clasificaciones["calcitriol"],

        "calcipotriol_clasificacion":
            clasificaciones["calcipotriol"]
    })

convergencia_df = pd.DataFrame(convergencia)

convergencia_df.to_csv(
    "CONTEXTO_genes_nucleo_convergencia.csv",
    index=False
)

# ================================================================
# 11. MOSTRAR RESULTADOS
# ================================================================

print("\n" + "=" * 90)
print("CONSISTENCIA POR TRATAMIENTO")
print("=" * 90)

for tratamiento in tratamientos_objetivo:

    sub = consistencia_df[
        consistencia_df["tratamiento"] == tratamiento
    ]

    if len(sub) == 0:
        continue

    print("\n" + tratamiento.upper())

    print(
        sub[
            [
                "gene",
                "n_contextos",
                "n_negativos",
                "n_positivos",
                "media_modz",
                "clasificacion"
            ]
        ]
        .sort_values(
            ["n_negativos", "media_modz"],
            ascending=[False, True]
        )
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

print("\n" + "=" * 90)
print("CONVERGENCIA ENTRE LOS TRES TRATAMIENTOS")
print("=" * 90)

if len(convergencia_df) > 0:

    print(
        convergencia_df
        .sort_values("gene")
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

else:
    print("No se pudo construir la tabla de convergencia.")

print("\nArchivos generados:")
print("- CONTEXTO_genes_nucleo_valores.csv")
print("- CONTEXTO_genes_nucleo_resumen.csv")
print("- CONTEXTO_genes_nucleo_signos.csv")
print("- CONTEXTO_genes_nucleo_consistencia.csv")
print("- CONTEXTO_genes_nucleo_convergencia.csv")

print("\nFIN")