import numpy as np
import pandas as pd

print("=" * 105)
print("HALLMARKS ROBUSTOS POR CONTEXTO - SENSIBILIDAD LANDMARK")
print("=" * 105)

# ================================================================
# 1. ARCHIVOS
# ================================================================

ARCHIVO_EXPR = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_GENES = "dataset_primario_genes.csv"
ARCHIVO_META = "dataset_primario_metadata.csv"
ARCHIVO_GMT = "data/h.all.v2026.1.Hs.symbols.gmt"
ARCHIVO_FULL = "HALLMARK_CONTEXTO_scores.csv"

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
    raise ValueError("No encontré pr_gene_symbol.")

if "pr_is_lm" not in genes.columns:
    raise ValueError("No encontré pr_is_lm.")

if "cell_id" not in meta.columns:
    raise ValueError("No encontré cell_id.")

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
# 4. IDENTIFICAR GENES LANDMARK
# ================================================================

lm_mask = (
    pd.to_numeric(
        genes["pr_is_lm"],
        errors="coerce"
    )
    .fillna(0)
    .astype(int)
    == 1
)

genes_lm = genes[lm_mask].copy()

print("\nGenes Landmark totales:", len(genes_lm))

# ================================================================
# 5. LEER HALLMARKS
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

        if nombre in HALLMARKS_OBJETIVO:
            hallmarks[nombre] = partes[2:]


if len(hallmarks) != len(HALLMARKS_OBJETIVO):
    raise ValueError(
        "No se encontraron los tres Hallmarks objetivo."
    )

# ================================================================
# 6. MAPA GEN -> ÍNDICE
# ================================================================

simbolos = genes["pr_gene_symbol"].astype(str).values

mapa_indices = {}

for i, simbolo in enumerate(simbolos):

    if simbolo not in mapa_indices:
        mapa_indices[simbolo] = i


landmark_symbols = set(
    genes_lm["pr_gene_symbol"].astype(str)
)

# ================================================================
# 7. COBERTURA LANDMARK
# ================================================================

cobertura = []
indices_hallmark_lm = {}

for hallmark in HALLMARKS_OBJETIVO:

    lista_genes = hallmarks[hallmark]

    presentes_lm = [
        g for g in lista_genes
        if g in mapa_indices
        and g in landmark_symbols
    ]

    indices = [
        mapa_indices[g]
        for g in presentes_lm
    ]

    indices_hallmark_lm[hallmark] = indices

    cobertura.append({
        "hallmark": hallmark,
        "genes_gmt": len(lista_genes),
        "genes_landmark_presentes": len(presentes_lm),
        "fraccion_landmark_sobre_gmt": (
            len(presentes_lm) / len(lista_genes)
            if len(lista_genes) > 0
            else np.nan
        )
    })


cobertura_df = pd.DataFrame(cobertura)

cobertura_df.to_csv(
    "HALLMARK_CONTEXTO_LM_cobertura.csv",
    index=False
)

print("\n" + "=" * 105)
print("COBERTURA LANDMARK")
print("=" * 105)

print(
    cobertura_df.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 8. CALCULAR SCORES LANDMARK POR FIRMA
# ================================================================

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

        idx = indices_hallmark_lm[hallmark]

        if len(idx) == 0:
            continue

        valores = X[i, idx]

        filas.append({
            "sig_id": sig_id,
            "cell_id": cell,
            "tratamiento": trat,
            "hallmark": hallmark,
            "n_genes_landmark": len(idx),
            "score_medio_LM": float(np.mean(valores)),
            "score_mediano_LM": float(np.median(valores)),
            "n_negativos_LM": int(np.sum(valores < 0)),
            "n_positivos_LM": int(np.sum(valores > 0)),
            "fraccion_negativa_LM": float(
                np.mean(valores < 0)
            )
        })


scores_lm = pd.DataFrame(filas)

scores_lm.to_csv(
    "HALLMARK_CONTEXTO_LM_scores.csv",
    index=False
)

# ================================================================
# 9. CONSENSO LANDMARK ENTRE LOS 3 TRATAMIENTOS
# ================================================================

consenso_lm = []

for cell in sorted(scores_lm["cell_id"].unique()):

    for hallmark in HALLMARKS_OBJETIVO:

        sub = scores_lm[
            (scores_lm["cell_id"] == cell)
            &
            (scores_lm["hallmark"] == hallmark)
        ]

        mapa = dict(
            zip(
                sub["tratamiento"],
                sub["score_medio_LM"]
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

        consenso_lm.append({
            "cell_id": cell,
            "hallmark": hallmark,
            "estradiol_LM": mapa["estradiol"],
            "calcitriol_LM": mapa["calcitriol"],
            "calcipotriol_LM": mapa["calcipotriol"],
            "patron_LM": patron
        })


consenso_lm_df = pd.DataFrame(consenso_lm)

consenso_lm_df.to_csv(
    "HALLMARK_CONTEXTO_LM_consenso.csv",
    index=False
)

# ================================================================
# 10. CARGAR RESULTADOS FULL DEL PASO 42
# ================================================================

full = pd.read_csv(ARCHIVO_FULL)

columnas_full = [
    "cell_id",
    "tratamiento",
    "hallmark",
    "score_medio"
]

faltan = [
    c for c in columnas_full
    if c not in full.columns
]

if faltan:
    raise ValueError(
        f"Faltan columnas en {ARCHIVO_FULL}: {faltan}"
    )

comparacion = full[
    columnas_full
].merge(
    scores_lm[
        [
            "cell_id",
            "tratamiento",
            "hallmark",
            "score_medio_LM",
            "fraccion_negativa_LM"
        ]
    ],
    on=[
        "cell_id",
        "tratamiento",
        "hallmark"
    ],
    how="inner"
)

# ================================================================
# 11. COMPARAR DIRECCIÓN FULL VS LANDMARK
# ================================================================

def signo(valor):

    if valor < 0:
        return "negativo"

    if valor > 0:
        return "positivo"

    return "cero"


comparacion["signo_Full"] = (
    comparacion["score_medio"].apply(signo)
)

comparacion["signo_LM"] = (
    comparacion["score_medio_LM"].apply(signo)
)

comparacion["misma_direccion"] = (
    comparacion["signo_Full"]
    ==
    comparacion["signo_LM"]
)

comparacion["diferencia_LM_menos_Full"] = (
    comparacion["score_medio_LM"]
    -
    comparacion["score_medio"]
)

comparacion.to_csv(
    "HALLMARK_CONTEXTO_LM_vs_FULL.csv",
    index=False
)

# ================================================================
# 12. RESUMEN DE CONCORDANCIA FULL-LM
# ================================================================

resumen_concordancia = (
    comparacion
    .groupby(
        "cell_id",
        as_index=False
    )
    .agg(
        comparaciones=("misma_direccion", "count"),
        misma_direccion=(
            "misma_direccion",
            "sum"
        )
    )
)

resumen_concordancia[
    "fraccion_misma_direccion"
] = (
    resumen_concordancia["misma_direccion"]
    /
    resumen_concordancia["comparaciones"]
)

resumen_concordancia.to_csv(
    "HALLMARK_CONTEXTO_LM_vs_FULL_resumen.csv",
    index=False
)

# ================================================================
# 13. RESUMEN DE CONVERGENCIA LANDMARK POR CONTEXTO
# ================================================================

resumen_lm = (
    consenso_lm_df
    .groupby(
        "cell_id",
        as_index=False
    )
    .agg(
        hallmarks_evaluados=("hallmark", "count"),
        hallmarks_negativos_3de3=(
            "patron_LM",
            lambda x: int(
                (x == "negativo_3de3").sum()
            )
        ),
        hallmarks_discordantes=(
            "patron_LM",
            lambda x: int(
                (x == "discordante").sum()
            )
        )
    )
)

resumen_lm[
    "fraccion_hallmarks_negativos_3de3"
] = (
    resumen_lm[
        "hallmarks_negativos_3de3"
    ]
    /
    resumen_lm[
        "hallmarks_evaluados"
    ]
)

resumen_lm.to_csv(
    "HALLMARK_CONTEXTO_LM_resumen.csv",
    index=False
)

# ================================================================
# 14. MOSTRAR SCORES LANDMARK
# ================================================================

print("\n" + "=" * 105)
print("SCORES LANDMARK POR CONTEXTO, TRATAMIENTO Y HALLMARK")
print("=" * 105)

for hallmark in HALLMARKS_OBJETIVO:

    print("\n" + hallmark)

    sub = scores_lm[
        scores_lm["hallmark"] == hallmark
    ]

    print(
        sub[
            [
                "cell_id",
                "tratamiento",
                "score_medio_LM",
                "score_mediano_LM",
                "fraccion_negativa_LM"
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
# 15. MOSTRAR CONVERGENCIA LANDMARK
# ================================================================

print("\n" + "=" * 105)
print("CONVERGENCIA LANDMARK ENTRE LOS TRES TRATAMIENTOS")
print("=" * 105)

print(
    consenso_lm_df
    .sort_values(
        ["cell_id", "hallmark"]
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 16. MOSTRAR FULL VS LANDMARK
# ================================================================

print("\n" + "=" * 105)
print("COMPARACIÓN FULL VS LANDMARK")
print("=" * 105)

print(
    comparacion[
        [
            "cell_id",
            "tratamiento",
            "hallmark",
            "score_medio",
            "score_medio_LM",
            "signo_Full",
            "signo_LM",
            "misma_direccion"
        ]
    ]
    .sort_values(
        [
            "cell_id",
            "hallmark",
            "tratamiento"
        ]
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 17. RESUMEN FINAL
# ================================================================

print("\n" + "=" * 105)
print("RESUMEN DE CONVERGENCIA LANDMARK POR CONTEXTO")
print("=" * 105)

print(
    resumen_lm
    .sort_values(
        "fraccion_hallmarks_negativos_3de3",
        ascending=False
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

print("\n" + "=" * 105)
print("CONCORDANCIA DIRECCIONAL FULL VS LANDMARK")
print("=" * 105)

print(
    resumen_concordancia
    .sort_values(
        "fraccion_misma_direccion",
        ascending=False
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

# ================================================================
# 18. NOTA
# ================================================================

print("\n" + "=" * 105)
print("NOTA METODOLÓGICA")
print("=" * 105)

print(
    "Este análisis utiliza exclusivamente genes Landmark directamente "
    "medidos por L1000."
)

print(
    "Los scores son medias descriptivas de MODZ dentro de cada Hallmark "
    "y no equivalen a NES ni a pruebas GSEA por contexto."
)

print(
    "La comparación Full vs Landmark evalúa robustez direccional frente "
    "a la inclusión de genes inferidos."
)

print(
    "Contexto celular y bloque experimental permanecen confundidos en "
    "el diseño primario; no se realizan pruebas inferenciales de línea celular."
)

print("\nArchivos generados:")
print("- HALLMARK_CONTEXTO_LM_cobertura.csv")
print("- HALLMARK_CONTEXTO_LM_scores.csv")
print("- HALLMARK_CONTEXTO_LM_consenso.csv")
print("- HALLMARK_CONTEXTO_LM_vs_FULL.csv")
print("- HALLMARK_CONTEXTO_LM_vs_FULL_resumen.csv")
print("- HALLMARK_CONTEXTO_LM_resumen.csv")

print("\nFIN")