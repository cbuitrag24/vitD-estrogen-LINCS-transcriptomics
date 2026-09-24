import numpy as np
import pandas as pd

print("=" * 95)
print("HETEROGENEIDAD CONTEXTUAL DEL NÚCLEO TRANSCRIPCIONAL COMPARTIDO")
print("=" * 95)

# ================================================================
# 1. ARCHIVOS
# ================================================================

ARCHIVO_EXPR = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_GENES = "dataset_primario_genes.csv"
ARCHIVO_META = "dataset_primario_metadata.csv"

# Genes Landmark recurrentes en >=2 Hallmarks robustos
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

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

# ================================================================
# 2. CARGAR DATOS
# ================================================================

npz = np.load(ARCHIVO_EXPR)
X = npz["expresion"]

genes = pd.read_csv(ARCHIVO_GENES)
meta = pd.read_csv(ARCHIVO_META)

print("\nDimensión expresión:", X.shape)
print("Dimensión metadata:", meta.shape)
print("Dimensión genes:", genes.shape)

if X.shape[0] != len(meta):
    raise ValueError(
        f"No coinciden las firmas: matriz={X.shape[0]}, metadata={len(meta)}"
    )

if X.shape[1] != len(genes):
    raise ValueError(
        f"No coinciden los genes: matriz={X.shape[1]}, genes={len(genes)}"
    )

# ================================================================
# 3. IDENTIFICAR TRATAMIENTO
# ================================================================

if "pert_iname" in meta.columns:
    col_trat = "pert_iname"
elif "tratamiento" in meta.columns:
    col_trat = "tratamiento"
else:
    raise ValueError(
        "No encontré pert_iname ni tratamiento en metadata."
    )

if "cell_id" not in meta.columns:
    raise ValueError("No encontré cell_id en metadata.")

if "pr_gene_symbol" not in genes.columns:
    raise ValueError("No encontré pr_gene_symbol.")

def normalizar_tratamiento(x):

    x = str(x).strip().lower()

    if "calcipotriol" in x:
        return "calcipotriol"

    if "calcitriol" in x:
        return "calcitriol"

    if "estradiol" in x:
        return "estradiol"

    return x

meta = meta.copy()

meta["tratamiento_normalizado"] = (
    meta[col_trat].apply(normalizar_tratamiento)
)

print(
    "\nContextos:",
    sorted(meta["cell_id"].astype(str).unique())
)

print(
    "Tratamientos:",
    sorted(meta["tratamiento_normalizado"].unique())
)

# ================================================================
# 4. LOCALIZAR GENES
# ================================================================

simbolos = genes["pr_gene_symbol"].astype(str).values

indices = {}
faltantes = []

for gene in GENES_NUCLEO:

    pos = np.where(simbolos == gene)[0]

    if len(pos) == 0:
        faltantes.append(gene)
    else:
        indices[gene] = pos[0]

print("\nGenes núcleo encontrados:", len(indices))
print(list(indices.keys()))

if faltantes:
    print("Genes faltantes:", faltantes)

if len(indices) != len(GENES_NUCLEO):
    print(
        "\nADVERTENCIA: no están presentes todos los genes núcleo."
    )

# ================================================================
# 5. CONSTRUIR TABLA LARGA
# ================================================================

filas = []

for i in range(len(meta)):

    cell = str(meta.iloc[i]["cell_id"])
    trat = meta.iloc[i]["tratamiento_normalizado"]

    for gene, j in indices.items():

        filas.append({
            "cell_id": cell,
            "tratamiento": trat,
            "gene": gene,
            "modz": float(X[i, j])
        })

datos = pd.DataFrame(filas)

datos.to_csv(
    "HETEROGENEIDAD_genes_nucleo_valores.csv",
    index=False
)

# ================================================================
# 6. RESUMEN POR CONTEXTO Y TRATAMIENTO
# ================================================================

resumen_contexto = (
    datos
    .groupby(
        ["cell_id", "tratamiento"],
        as_index=False
    )
    .agg(
        media_nucleo=("modz", "mean"),
        mediana_nucleo=("modz", "median"),
        minimo=("modz", "min"),
        maximo=("modz", "max"),
        n_genes=("gene", "nunique"),
        n_negativos=("modz", lambda x: int((x < 0).sum())),
        n_positivos=("modz", lambda x: int((x > 0).sum()))
    )
)

resumen_contexto["fraccion_negativa"] = (
    resumen_contexto["n_negativos"]
    / resumen_contexto["n_genes"]
)

resumen_contexto["fraccion_positiva"] = (
    resumen_contexto["n_positivos"]
    / resumen_contexto["n_genes"]
)

resumen_contexto.to_csv(
    "HETEROGENEIDAD_resumen_contexto_tratamiento.csv",
    index=False
)

# ================================================================
# 7. CONSENSO DE LOS 3 TRATAMIENTOS DENTRO DE CADA CONTEXTO
# ================================================================

consenso = []

for cell in sorted(datos["cell_id"].unique()):

    sub_cell = datos[
        datos["cell_id"] == cell
    ]

    for gene in indices:

        sub_gene = sub_cell[
            sub_cell["gene"] == gene
        ]

        mapa = dict(
            zip(
                sub_gene["tratamiento"],
                sub_gene["modz"]
            )
        )

        if not all(
            t in mapa
            for t in TRATAMIENTOS
        ):
            continue

        valores = [
            mapa[t]
            for t in TRATAMIENTOS
        ]

        n_neg = sum(v < 0 for v in valores)
        n_pos = sum(v > 0 for v in valores)

        if n_neg == 3:
            patron = "negativo_3de3"
        elif n_pos == 3:
            patron = "positivo_3de3"
        else:
            patron = "discordante"

        consenso.append({
            "cell_id": cell,
            "gene": gene,
            "estradiol": mapa["estradiol"],
            "calcitriol": mapa["calcitriol"],
            "calcipotriol": mapa["calcipotriol"],
            "n_tratamientos_negativos": n_neg,
            "n_tratamientos_positivos": n_pos,
            "patron": patron
        })

consenso_df = pd.DataFrame(consenso)

consenso_df.to_csv(
    "HETEROGENEIDAD_consenso_tres_tratamientos_por_contexto.csv",
    index=False
)

# ================================================================
# 8. RESUMEN DE CONCORDANCIA POR CONTEXTO
# ================================================================

resumen_consenso = (
    consenso_df
    .groupby(
        "cell_id",
        as_index=False
    )
    .agg(
        genes_evaluados=("gene", "nunique"),
        genes_negativos_3de3=(
            "patron",
            lambda x: int((x == "negativo_3de3").sum())
        ),
        genes_positivos_3de3=(
            "patron",
            lambda x: int((x == "positivo_3de3").sum())
        ),
        genes_discordantes=(
            "patron",
            lambda x: int((x == "discordante").sum())
        )
    )
)

resumen_consenso["fraccion_negativos_3de3"] = (
    resumen_consenso["genes_negativos_3de3"]
    / resumen_consenso["genes_evaluados"]
)

resumen_consenso.to_csv(
    "HETEROGENEIDAD_resumen_concordancia_contextos.csv",
    index=False
)

# ================================================================
# 9. IDENTIFICAR EXCEPCIONES POSITIVAS
# ================================================================

excepciones = datos[
    datos["modz"] > 0
].copy()

excepciones = excepciones.sort_values(
    ["cell_id", "tratamiento", "modz"],
    ascending=[True, True, False]
)

excepciones.to_csv(
    "HETEROGENEIDAD_excepciones_positivas.csv",
    index=False
)

# ================================================================
# 10. CONTAR EXCEPCIONES POR CONTEXTO
# ================================================================

conteo_excepciones = (
    excepciones
    .groupby(
        ["cell_id", "tratamiento"],
        as_index=False
    )
    .agg(
        n_excepciones_positivas=("gene", "count")
    )
)

base = (
    datos[
        ["cell_id", "tratamiento"]
    ]
    .drop_duplicates()
)

conteo_excepciones = base.merge(
    conteo_excepciones,
    on=["cell_id", "tratamiento"],
    how="left"
)

conteo_excepciones[
    "n_excepciones_positivas"
] = (
    conteo_excepciones[
        "n_excepciones_positivas"
    ]
    .fillna(0)
    .astype(int)
)

conteo_excepciones.to_csv(
    "HETEROGENEIDAD_conteo_excepciones.csv",
    index=False
)

# ================================================================
# 11. AMPLITUD DE HETEROGENEIDAD POR GEN
# ================================================================

heterogeneidad_gen = (
    datos
    .groupby(
        ["tratamiento", "gene"],
        as_index=False
    )
    .agg(
        minimo_modz=("modz", "min"),
        maximo_modz=("modz", "max"),
        media_modz=("modz", "mean"),
        desviacion_contextual=("modz", "std")
    )
)

heterogeneidad_gen["rango_contextual"] = (
    heterogeneidad_gen["maximo_modz"]
    - heterogeneidad_gen["minimo_modz"]
)

heterogeneidad_gen = (
    heterogeneidad_gen
    .sort_values(
        ["tratamiento", "rango_contextual"],
        ascending=[True, False]
    )
)

heterogeneidad_gen.to_csv(
    "HETEROGENEIDAD_por_gen.csv",
    index=False
)

# ================================================================
# 12. MOSTRAR RESULTADOS
# ================================================================

print("\n" + "=" * 95)
print("RESUMEN DEL NÚCLEO POR CONTEXTO Y TRATAMIENTO")
print("=" * 95)

print(
    resumen_contexto[
        [
            "cell_id",
            "tratamiento",
            "media_nucleo",
            "mediana_nucleo",
            "n_negativos",
            "n_positivos",
            "fraccion_negativa"
        ]
    ]
    .sort_values(
        ["cell_id", "tratamiento"]
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

print("\n" + "=" * 95)
print("CONCORDANCIA DE LOS 3 TRATAMIENTOS DENTRO DE CADA CONTEXTO")
print("=" * 95)

print(
    resumen_consenso
    .sort_values(
        "fraccion_negativos_3de3",
        ascending=False
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

print("\n" + "=" * 95)
print("EXCEPCIONES POSITIVAS")
print("=" * 95)

if len(excepciones) > 0:

    print(
        excepciones[
            [
                "cell_id",
                "tratamiento",
                "gene",
                "modz"
            ]
        ]
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

else:
    print("No se encontraron excepciones positivas.")

print("\n" + "=" * 95)
print("GENES CON MAYOR HETEROGENEIDAD CONTEXTUAL")
print("=" * 95)

for tratamiento in TRATAMIENTOS:

    sub = heterogeneidad_gen[
        heterogeneidad_gen["tratamiento"] == tratamiento
    ].head(5)

    print("\n" + tratamiento.upper())

    print(
        sub[
            [
                "gene",
                "minimo_modz",
                "maximo_modz",
                "rango_contextual",
                "desviacion_contextual"
            ]
        ]
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

# ================================================================
# 13. NOTA METODOLÓGICA
# ================================================================

print("\n" + "=" * 95)
print("NOTA METODOLÓGICA")
print("=" * 95)

print(
    "Este análisis es descriptivo. Cada contexto celular del diseño "
    "primario está asociado a un bloque experimental específico."
)

print(
    "Por ello, contexto celular y bloque experimental están confundidos "
    "y no se interpreta este análisis como una prueba inferencial "
    "independiente del efecto de línea celular."
)

print("\nArchivos generados:")
print("- HETEROGENEIDAD_genes_nucleo_valores.csv")
print("- HETEROGENEIDAD_resumen_contexto_tratamiento.csv")
print("- HETEROGENEIDAD_consenso_tres_tratamientos_por_contexto.csv")
print("- HETEROGENEIDAD_resumen_concordancia_contextos.csv")
print("- HETEROGENEIDAD_excepciones_positivas.csv")
print("- HETEROGENEIDAD_conteo_excepciones.csv")
print("- HETEROGENEIDAD_por_gen.csv")

print("\nFIN")