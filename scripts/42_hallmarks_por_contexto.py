import numpy as np
import pandas as pd

print("=" * 100)
print("HALLMARKS ROBUSTOS POR CONTEXTO")
print("=" * 100)

# ================================================================
# 1. ARCHIVOS
# ================================================================

ARCHIVO_EXPR = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_GENES = "dataset_primario_genes.csv"
ARCHIVO_META = "dataset_primario_metadata.csv"
ARCHIVO_GMT = "data/h.all.v2026.1.Hs.symbols.gmt"

HALLMARKS_OBJETIVO = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1"
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
        f"No coinciden firmas: matriz={X.shape[0]}, metadata={len(meta)}"
    )

if X.shape[1] != len(genes):
    raise ValueError(
        f"No coinciden genes: matriz={X.shape[1]}, genes={len(genes)}"
    )

if "pr_gene_symbol" not in genes.columns:
    raise ValueError("No encontré pr_gene_symbol en genes.")

if "cell_id" not in meta.columns:
    raise ValueError("No encontré cell_id en metadata.")

if "pert_iname" in meta.columns:
    col_trat = "pert_iname"
elif "tratamiento" in meta.columns:
    col_trat = "tratamiento"
else:
    raise ValueError("No encontré columna de tratamiento.")

# ================================================================
# 3. NORMALIZAR TRATAMIENTOS
# ================================================================

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
# 4. LEER GMT
# ================================================================

hallmarks = {}

with open(
    ARCHIVO_GMT,
    "r",
    encoding="utf-8"
) as f:

    for linea in f:

        partes = linea.rstrip("\n").split("\t")

        if len(partes) < 3:
            continue

        nombre = partes[0]
        genes_set = partes[2:]

        if nombre in HALLMARKS_OBJETIVO:
            hallmarks[nombre] = genes_set


print("\nHallmarks encontrados:", len(hallmarks))

for h in HALLMARKS_OBJETIVO:

    if h in hallmarks:
        print(
            h,
            "->",
            len(hallmarks[h]),
            "genes en GMT"
        )
    else:
        print(
            "NO ENCONTRADO:",
            h
        )

if len(hallmarks) != len(HALLMARKS_OBJETIVO):
    raise ValueError(
        "No se encontraron los tres Hallmarks objetivo en el GMT."
    )

# ================================================================
# 5. MAPA DE GENES DEL DATASET
# ================================================================

simbolos = (
    genes["pr_gene_symbol"]
    .astype(str)
    .values
)

mapa_indices = {}

for i, simbolo in enumerate(simbolos):

    if simbolo not in mapa_indices:
        mapa_indices[simbolo] = i

# ================================================================
# 6. COBERTURA DE CADA HALLMARK
# ================================================================

cobertura = []
indices_hallmark = {}

for hallmark, lista_genes in hallmarks.items():

    presentes = [
        g for g in lista_genes
        if g in mapa_indices
    ]

    indices = [
        mapa_indices[g]
        for g in presentes
    ]

    indices_hallmark[hallmark] = indices

    cobertura.append({
        "hallmark": hallmark,
        "genes_gmt": len(lista_genes),
        "genes_presentes_L1000": len(presentes),
        "fraccion_cobertura": (
            len(presentes) / len(lista_genes)
            if len(lista_genes) > 0
            else np.nan
        )
    })


cobertura_df = pd.DataFrame(cobertura)

cobertura_df.to_csv(
    "HALLMARK_CONTEXTO_cobertura.csv",
    index=False
)

print("\n" + "=" * 100)
print("COBERTURA")
print("=" * 100)

print(
    cobertura_df.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 7. CALCULAR SCORE POR FIRMA
# ================================================================

#
# Score descriptivo:
# media de los MODZ de todos los genes del Hallmark
#
# Como los tres programas fueron negativos en el análisis GSEA,
# un score más negativo indica mayor supresión transcriptómica media.
#

filas = []

for i in range(len(meta)):

    cell = str(meta.iloc[i]["cell_id"])
    trat = meta.iloc[i]["tratamiento_normalizado"]

    sig_id = (
        meta.iloc[i]["sig_id"]
        if "sig_id" in meta.columns
        else f"firma_{i}"
    )

    for hallmark in HALLMARKS_OBJETIVO:

        idx = indices_hallmark[hallmark]

        valores = X[i, idx]

        filas.append({
            "sig_id": sig_id,
            "cell_id": cell,
            "tratamiento": trat,
            "hallmark": hallmark,
            "n_genes": len(idx),
            "score_medio": float(np.mean(valores)),
            "score_mediano": float(np.median(valores)),
            "n_negativos": int(np.sum(valores < 0)),
            "n_positivos": int(np.sum(valores > 0)),
            "fraccion_negativa": float(
                np.mean(valores < 0)
            )
        })


scores = pd.DataFrame(filas)

scores.to_csv(
    "HALLMARK_CONTEXTO_scores.csv",
    index=False
)

# ================================================================
# 8. CONSENSO ENTRE LOS TRES TRATAMIENTOS
# ================================================================

consenso = []

for cell in sorted(scores["cell_id"].unique()):

    for hallmark in HALLMARKS_OBJETIVO:

        sub = scores[
            (scores["cell_id"] == cell)
            &
            (scores["hallmark"] == hallmark)
        ]

        mapa = dict(
            zip(
                sub["tratamiento"],
                sub["score_medio"]
            )
        )

        if not all(
            t in mapa
            for t in TRATAMIENTOS
        ):
            continue

        vals = [
            mapa[t]
            for t in TRATAMIENTOS
        ]

        n_neg = sum(v < 0 for v in vals)
        n_pos = sum(v > 0 for v in vals)

        if n_neg == 3:
            patron = "negativo_3de3"
        elif n_pos == 3:
            patron = "positivo_3de3"
        else:
            patron = "discordante"

        consenso.append({
            "cell_id": cell,
            "hallmark": hallmark,
            "estradiol": mapa["estradiol"],
            "calcitriol": mapa["calcitriol"],
            "calcipotriol": mapa["calcipotriol"],
            "n_negativos": n_neg,
            "n_positivos": n_pos,
            "patron": patron
        })


consenso_df = pd.DataFrame(consenso)

consenso_df.to_csv(
    "HALLMARK_CONTEXTO_consenso.csv",
    index=False
)

# ================================================================
# 9. RESUMEN POR CONTEXTO
# ================================================================

resumen_contexto = (
    consenso_df
    .groupby(
        "cell_id",
        as_index=False
    )
    .agg(
        hallmarks_evaluados=("hallmark", "count"),
        hallmarks_negativos_3de3=(
            "patron",
            lambda x: int(
                (x == "negativo_3de3").sum()
            )
        ),
        hallmarks_discordantes=(
            "patron",
            lambda x: int(
                (x == "discordante").sum()
            )
        )
    )
)

resumen_contexto[
    "fraccion_hallmarks_negativos_3de3"
] = (
    resumen_contexto[
        "hallmarks_negativos_3de3"
    ]
    /
    resumen_contexto[
        "hallmarks_evaluados"
    ]
)

resumen_contexto.to_csv(
    "HALLMARK_CONTEXTO_resumen.csv",
    index=False
)

# ================================================================
# 10. MATRIZ PARA VISUALIZAR
# ================================================================

matriz = (
    scores
    .pivot_table(
        index=[
            "hallmark",
            "cell_id"
        ],
        columns="tratamiento",
        values="score_medio",
        aggfunc="mean"
    )
    .reset_index()
)

matriz.to_csv(
    "HALLMARK_CONTEXTO_matriz.csv",
    index=False
)

# ================================================================
# 11. MOSTRAR SCORES
# ================================================================

print("\n" + "=" * 100)
print("SCORES MEDIOS POR CONTEXTO, TRATAMIENTO Y HALLMARK")
print("=" * 100)

for hallmark in HALLMARKS_OBJETIVO:

    print("\n" + hallmark)

    sub = scores[
        scores["hallmark"] == hallmark
    ]

    print(
        sub[
            [
                "cell_id",
                "tratamiento",
                "score_medio",
                "score_mediano",
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

# ================================================================
# 12. MOSTRAR CONSENSO
# ================================================================

print("\n" + "=" * 100)
print("CONSENSO ENTRE LOS TRES TRATAMIENTOS")
print("=" * 100)

print(
    consenso_df[
        [
            "cell_id",
            "hallmark",
            "estradiol",
            "calcitriol",
            "calcipotriol",
            "patron"
        ]
    ]
    .sort_values(
        ["cell_id", "hallmark"]
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 13. RESUMEN FINAL POR CONTEXTO
# ================================================================

print("\n" + "=" * 100)
print("RESUMEN DE CONVERGENCIA DE HALLMARKS POR CONTEXTO")
print("=" * 100)

print(
    resumen_contexto
    .sort_values(
        "fraccion_hallmarks_negativos_3de3",
        ascending=False
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 14. NOTA METODOLÓGICA
# ================================================================

print("\n" + "=" * 100)
print("NOTA METODOLÓGICA")
print("=" * 100)

print(
    "Los scores por contexto son descriptivos y corresponden a la "
    "media de los valores MODZ de los genes presentes de cada Hallmark."
)

print(
    "No equivalen al NES del análisis GSEA global y no deben "
    "interpretarse como una prueba de enriquecimiento por contexto."
)

print(
    "Cada contexto celular está asociado a un único bloque experimental; "
    "por ello, contexto y bloque están confundidos en el diseño primario."
)

print("\nArchivos generados:")
print("- HALLMARK_CONTEXTO_cobertura.csv")
print("- HALLMARK_CONTEXTO_scores.csv")
print("- HALLMARK_CONTEXTO_consenso.csv")
print("- HALLMARK_CONTEXTO_resumen.csv")
print("- HALLMARK_CONTEXTO_matriz.csv")

print("\nFIN")