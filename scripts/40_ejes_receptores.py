import numpy as np
import pandas as pd

print("=" * 95)
print("ANÁLISIS DE EJES RECEPTORES: VITAMINA D Y ESTRÓGENO")
print("=" * 95)

# ================================================================
# 1. ARCHIVOS
# ================================================================

ARCHIVO_EXPR = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_GENES = "dataset_primario_genes.csv"
ARCHIVO_META = "dataset_primario_metadata.csv"

# ================================================================
# 2. GENES DE INTERÉS
# ================================================================

EJE_VITAMINA_D = [
    "VDR",
    "RXRA",
    "RXRB",
    "RXRG",
    "CYP24A1",
    "CYP27B1"
]

EJE_ESTROGENO = [
    "ESR1",
    "ESR2",
    "GPER1"
]

GENES_INTERES = EJE_VITAMINA_D + EJE_ESTROGENO

# ================================================================
# 3. CARGAR DATOS
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
        f"No coinciden genes: matriz={X.shape[1]}, tabla={len(genes)}"
    )

# ================================================================
# 4. COMPROBAR COLUMNAS
# ================================================================

if "pr_gene_symbol" not in genes.columns:
    raise ValueError("Falta pr_gene_symbol en la tabla de genes.")

if "cell_id" not in meta.columns:
    raise ValueError("Falta cell_id en metadata.")

if "pert_iname" in meta.columns:
    col_trat = "pert_iname"
elif "tratamiento" in meta.columns:
    col_trat = "tratamiento"
else:
    raise ValueError(
        "No encontré pert_iname ni tratamiento en metadata."
    )

# ================================================================
# 5. NORMALIZAR TRATAMIENTOS
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
    "\nTratamientos:",
    sorted(meta["tratamiento_normalizado"].unique())
)

print(
    "Contextos:",
    sorted(meta["cell_id"].astype(str).unique())
)

# ================================================================
# 6. LOCALIZAR GENES
# ================================================================

simbolos = genes["pr_gene_symbol"].astype(str).values

indices = {}
faltantes = []

for gene in GENES_INTERES:

    encontrados = np.where(simbolos == gene)[0]

    if len(encontrados) == 0:
        faltantes.append(gene)
    else:
        indices[gene] = encontrados[0]


print("\n" + "=" * 95)
print("COBERTURA DE GENES")
print("=" * 95)

print("\nGenes encontrados:")
print(list(indices.keys()))

if faltantes:
    print("\nGenes no encontrados:")
    print(faltantes)

# ================================================================
# 7. INFORMACIÓN LANDMARK / INFERIDO
# ================================================================

cobertura = []

for gene in GENES_INTERES:

    if gene not in indices:

        cobertura.append({
            "gene": gene,
            "eje": (
                "vitamina_D"
                if gene in EJE_VITAMINA_D
                else "estrogeno"
            ),
            "presente": False,
            "pr_is_lm": np.nan,
            "tipo_medicion": "no_presente"
        })

        continue

    j = indices[gene]

    if "pr_is_lm" in genes.columns:

        valor_lm = genes.iloc[j]["pr_is_lm"]

        try:
            es_landmark = int(valor_lm) == 1
        except:
            es_landmark = False

    else:

        valor_lm = np.nan
        es_landmark = False

    cobertura.append({
        "gene": gene,
        "eje": (
            "vitamina_D"
            if gene in EJE_VITAMINA_D
            else "estrogeno"
        ),
        "presente": True,
        "pr_is_lm": valor_lm,
        "tipo_medicion": (
            "landmark_directo"
            if es_landmark
            else "inferido"
        )
    })


cobertura_df = pd.DataFrame(cobertura)

cobertura_df.to_csv(
    "RECEPTORES_cobertura_genes.csv",
    index=False
)

print("\nCobertura:")
print(
    cobertura_df.to_string(
        index=False
    )
)

# ================================================================
# 8. EXTRAER MODZ POR FIRMA
# ================================================================

filas = []

for i in range(len(meta)):

    contexto = meta.iloc[i]["cell_id"]
    tratamiento = meta.iloc[i]["tratamiento_normalizado"]

    sig_id = (
        meta.iloc[i]["sig_id"]
        if "sig_id" in meta.columns
        else f"firma_{i}"
    )

    for gene, j in indices.items():

        eje = (
            "vitamina_D"
            if gene in EJE_VITAMINA_D
            else "estrogeno"
        )

        filas.append({
            "sig_id": sig_id,
            "cell_id": contexto,
            "tratamiento": tratamiento,
            "eje": eje,
            "gene": gene,
            "modz": float(X[i, j])
        })


valores = pd.DataFrame(filas)

valores.to_csv(
    "RECEPTORES_valores_por_contexto.csv",
    index=False
)

# ================================================================
# 9. MATRIZ CONTEXTO × TRATAMIENTO
# ================================================================

matriz = (
    valores
    .pivot_table(
        index=["gene", "eje"],
        columns=["tratamiento", "cell_id"],
        values="modz",
        aggfunc="mean"
    )
)

matriz.to_csv(
    "RECEPTORES_matriz_contextos.csv"
)

# ================================================================
# 10. RESUMEN POR TRATAMIENTO Y GEN
# ================================================================

resumen = (
    valores
    .groupby(
        ["tratamiento", "eje", "gene"],
        as_index=False
    )
    .agg(
        media_modz=("modz", "mean"),
        mediana_modz=("modz", "median"),
        minimo=("modz", "min"),
        maximo=("modz", "max"),
        n_contextos=("cell_id", "nunique")
    )
)

# Conteo de signos
signos = (
    valores
    .assign(
        negativo=lambda x: x["modz"] < 0,
        positivo=lambda x: x["modz"] > 0
    )
    .groupby(
        ["tratamiento", "eje", "gene"],
        as_index=False
    )
    .agg(
        n_negativos=("negativo", "sum"),
        n_positivos=("positivo", "sum")
    )
)

resumen = resumen.merge(
    signos,
    on=["tratamiento", "eje", "gene"],
    how="left"
)

# ================================================================
# 11. CLASIFICACIÓN DE CONSISTENCIA
# ================================================================

def clasificar(row):

    n = int(row["n_contextos"])
    neg = int(row["n_negativos"])
    pos = int(row["n_positivos"])

    if neg == n:
        return f"negativo_{n}de{n}"

    if pos == n:
        return f"positivo_{n}de{n}"

    if n == 4 and neg == 3:
        return "predominio_negativo_3de4"

    if n == 4 and pos == 3:
        return "predominio_positivo_3de4"

    return "heterogeneo"


resumen["clasificacion"] = (
    resumen.apply(
        clasificar,
        axis=1
    )
)

resumen.to_csv(
    "RECEPTORES_resumen_tratamientos.csv",
    index=False
)

# ================================================================
# 12. RESPUESTA DE CADA GEN EN CADA CONTEXTO
# ================================================================

contexto = (
    valores
    .pivot_table(
        index=["gene", "eje", "cell_id"],
        columns="tratamiento",
        values="modz",
        aggfunc="mean"
    )
    .reset_index()
)

contexto.to_csv(
    "RECEPTORES_comparacion_por_contexto.csv",
    index=False
)

# ================================================================
# 13. DIFERENCIAS ENTRE ESTRADIOL Y VITAMINA D
# ================================================================

comparaciones = []

for gene in indices:

    sub = contexto[
        contexto["gene"] == gene
    ]

    for _, row in sub.iterrows():

        e2 = row.get("estradiol", np.nan)
        calc = row.get("calcitriol", np.nan)
        calcip = row.get("calcipotriol", np.nan)

        promedio_vitd = np.nanmean(
            [calc, calcip]
        )

        comparaciones.append({
            "gene": gene,
            "eje": row["eje"],
            "cell_id": row["cell_id"],
            "estradiol_modz": e2,
            "calcitriol_modz": calc,
            "calcipotriol_modz": calcip,
            "promedio_vitaminaD": promedio_vitd,
            "E2_menos_VitD": (
                e2 - promedio_vitd
                if not np.isnan(e2)
                else np.nan
            )
        })


comparaciones_df = pd.DataFrame(comparaciones)

comparaciones_df.to_csv(
    "RECEPTORES_E2_vs_VitD_por_contexto.csv",
    index=False
)

# ================================================================
# 14. MOSTRAR RESULTADOS IMPORTANTES
# ================================================================

print("\n" + "=" * 95)
print("RESPUESTA TRANSCRIPTÓMICA DE LOS GENES DEL EJE RECEPTOR")
print("=" * 95)

for tratamiento in [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]:

    sub = resumen[
        resumen["tratamiento"] == tratamiento
    ]

    if len(sub) == 0:
        continue

    print("\n" + tratamiento.upper())

    print(
        sub[
            [
                "eje",
                "gene",
                "media_modz",
                "n_negativos",
                "n_positivos",
                "clasificacion"
            ]
        ]
        .sort_values(
            ["eje", "gene"]
        )
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

# ================================================================
# 15. MOSTRAR CONTEXTOS INDIVIDUALES
# ================================================================

print("\n" + "=" * 95)
print("VALORES POR CONTEXTO")
print("=" * 95)

for gene in GENES_INTERES:

    if gene not in indices:
        continue

    print(f"\n{gene}")

    sub = contexto[
        contexto["gene"] == gene
    ]

    columnas = [
        c for c in [
            "cell_id",
            "estradiol",
            "calcitriol",
            "calcipotriol"
        ]
        if c in sub.columns
    ]

    print(
        sub[columnas]
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

# ================================================================
# 16. ADVERTENCIA METODOLÓGICA
# ================================================================

print("\n" + "=" * 95)
print("NOTA METODOLÓGICA")
print("=" * 95)

print(
    "Los valores analizados son firmas diferenciales MODZ de LINCS L1000."
)

print(
    "NO representan expresión basal del receptor y NO constituyen "
    "una medida directa de actividad funcional del receptor."
)

print(
    "Los genes clasificados como no-landmark corresponden a valores "
    "inferidos en el espacio L1000."
)

print("\nArchivos generados:")
print("- RECEPTORES_cobertura_genes.csv")
print("- RECEPTORES_valores_por_contexto.csv")
print("- RECEPTORES_matriz_contextos.csv")
print("- RECEPTORES_resumen_tratamientos.csv")
print("- RECEPTORES_comparacion_por_contexto.csv")
print("- RECEPTORES_E2_vs_VitD_por_contexto.csv")

print("\nFIN")