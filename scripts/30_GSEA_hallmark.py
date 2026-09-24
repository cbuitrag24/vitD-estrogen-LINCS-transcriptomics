import numpy as np
import pandas as pd
from pathlib import Path

# ================================================================
# CONFIGURACIÓN
# ================================================================

N_PERM = 5000
SEED = 20260903
MIN_GENES = 15
WEIGHT = 1.0

TRATAMIENTOS = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

GMT_FILE = "h.all.v2026.1.Hs.symbols.gmt"

print("=" * 80)
print("GSEA PRERANKED - HALLMARK")
print("=" * 80)

print("\nConfiguración:")
print("Permutaciones por Hallmark:", N_PERM)
print("Semilla aleatoria:", SEED)
print("Mínimo de genes:", MIN_GENES)
print("Exponente de ponderación:", WEIGHT)


# ================================================================
# FUNCIONES
# ================================================================

def calcular_es(scores, hit_indices, weight=1.0):
    """
    Calcula el enrichment score (ES) de un gene set
    sobre un ranking ordenado de mayor a menor.
    """

    n_total = len(scores)
    hit_indices = np.asarray(hit_indices, dtype=int)
    n_hits = len(hit_indices)

    if n_hits == 0 or n_hits == n_total:
        return np.nan

    hit_mask = np.zeros(n_total, dtype=bool)
    hit_mask[hit_indices] = True

    pesos_hits = np.abs(scores[hit_indices]) ** weight
    suma_pesos = pesos_hits.sum()

    if suma_pesos == 0:
        pesos_hits = np.ones(n_hits, dtype=float)
        suma_pesos = float(n_hits)

    running = np.empty(n_total, dtype=float)

    incremento = np.zeros(n_total, dtype=float)
    incremento[hit_indices] = pesos_hits / suma_pesos

    decremento = np.where(
        hit_mask,
        0.0,
        1.0 / (n_total - n_hits)
    )

    running = np.cumsum(incremento - decremento)

    max_es = running.max()
    min_es = running.min()

    if abs(max_es) >= abs(min_es):
        return max_es
    else:
        return min_es


def bh_adjust(pvalues):
    """
    Corrección Benjamini-Hochberg.
    """
    pvalues = np.asarray(pvalues, dtype=float)

    n = len(pvalues)
    orden = np.argsort(pvalues)
    p_ordenados = pvalues[orden]

    q = np.empty(n, dtype=float)

    valores = p_ordenados * n / np.arange(1, n + 1)
    valores = np.minimum.accumulate(valores[::-1])[::-1]
    valores = np.clip(valores, 0, 1)

    q[orden] = valores

    return q


def leer_gmt(archivo):
    hallmarks = {}

    with open(archivo, "r", encoding="utf-8") as f:
        for linea in f:
            partes = linea.rstrip("\n").split("\t")

            if len(partes) < 3:
                continue

            nombre = partes[0]
            genes = set(partes[2:])

            hallmarks[nombre] = genes

    return hallmarks


# ================================================================
# LEER HALLMARK
# ================================================================

if not Path(GMT_FILE).exists():
    raise FileNotFoundError(
        f"No se encontró el archivo {GMT_FILE}"
    )

hallmarks = leer_gmt(GMT_FILE)

print("\nHallmark gene sets leídos:", len(hallmarks))

if len(hallmarks) != 50:
    raise ValueError(
        f"Se esperaban 50 Hallmarks y se encontraron {len(hallmarks)}"
    )


# ================================================================
# GENERADOR ALEATORIO
# ================================================================

rng = np.random.default_rng(SEED)

todos_resultados = []


# ================================================================
# ANALIZAR CADA TRATAMIENTO
# ================================================================

for tratamiento in TRATAMIENTOS:

    print("\n" + "=" * 80)
    print("TRATAMIENTO:", tratamiento.upper())
    print("=" * 80)

    archivo_ranking = f"RANKING_{tratamiento}.csv"

    if not Path(archivo_ranking).exists():
        raise FileNotFoundError(
            f"No se encontró {archivo_ranking}"
        )

    ranking = pd.read_csv(archivo_ranking)

    columnas_necesarias = {
        "pr_gene_symbol",
        "mean_modz"
    }

    if not columnas_necesarias.issubset(ranking.columns):
        raise ValueError(
            f"{archivo_ranking} no contiene las columnas esperadas"
        )

    ranking = ranking[
        ["pr_gene_symbol", "mean_modz"]
    ].copy()

    ranking["pr_gene_symbol"] = (
        ranking["pr_gene_symbol"]
        .astype(str)
        .str.strip()
    )

    ranking["mean_modz"] = pd.to_numeric(
        ranking["mean_modz"],
        errors="coerce"
    )

    ranking = ranking.dropna()

    # Comprobar duplicados
    duplicados = ranking["pr_gene_symbol"].duplicated().sum()

    if duplicados > 0:
        raise ValueError(
            f"Se encontraron {duplicados} símbolos génicos duplicados "
            f"en {archivo_ranking}. Se detiene el análisis."
        )

    # Orden reproducible:
    # primero mean_modz descendente,
    # luego símbolo génico alfabético para resolver empates
    ranking = ranking.sort_values(
        ["mean_modz", "pr_gene_symbol"],
        ascending=[False, True],
        kind="mergesort"
    ).reset_index(drop=True)

    genes_ranking = ranking["pr_gene_symbol"].to_numpy()
    scores = ranking["mean_modz"].to_numpy(dtype=float)

    n_total = len(ranking)

    print("Genes en ranking:", n_total)

    mapa_posicion = {
        gen: i
        for i, gen in enumerate(genes_ranking)
    }

    resultados_tratamiento = []

    for contador, (hallmark, genes_set) in enumerate(
        hallmarks.items(),
        start=1
    ):

        genes_presentes = [
            gen
            for gen in genes_set
            if gen in mapa_posicion
        ]

        n_hits = len(genes_presentes)

        if n_hits < MIN_GENES:
            continue

        hit_indices = np.array(
            [mapa_posicion[gen] for gen in genes_presentes],
            dtype=int
        )

        es_obs = calcular_es(
            scores,
            hit_indices,
            weight=WEIGHT
        )

        # --------------------------------------------------------
        # Distribución nula:
        # conjuntos aleatorios del mismo tamaño
        # --------------------------------------------------------

        null_es = np.empty(N_PERM, dtype=float)

        for i in range(N_PERM):

            indices_perm = rng.choice(
                n_total,
                size=n_hits,
                replace=False
            )

            null_es[i] = calcular_es(
                scores,
                indices_perm,
                weight=WEIGHT
            )

        # --------------------------------------------------------
        # NES y p empírico, separados por signo
        # --------------------------------------------------------

        if es_obs >= 0:

            null_mismo_signo = null_es[null_es >= 0]

            if len(null_mismo_signo) == 0:
                nes = np.nan
            else:
                denominador = np.mean(null_mismo_signo)
                nes = es_obs / denominador

            p_empirico = (
                np.sum(null_es >= es_obs) + 1
            ) / (N_PERM + 1)

        else:

            null_mismo_signo = null_es[null_es < 0]

            if len(null_mismo_signo) == 0:
                nes = np.nan
            else:
                denominador = np.mean(np.abs(null_mismo_signo))
                nes = es_obs / denominador

            p_empirico = (
                np.sum(null_es <= es_obs) + 1
            ) / (N_PERM + 1)

        resultados_tratamiento.append({
            "tratamiento": tratamiento,
            "hallmark": hallmark,
            "n_genes_hallmark_en_ranking": n_hits,
            "ES": es_obs,
            "NES": nes,
            "p_empirico": p_empirico
        })

        print(
            f"{contador:02d}/50  "
            f"{hallmark:<45} "
            f"n={n_hits:3d}  "
            f"ES={es_obs: .4f}  "
            f"NES={nes: .4f}  "
            f"p={p_empirico:.6f}"
        )

    resultados_tratamiento = pd.DataFrame(
        resultados_tratamiento
    )

    # BH dentro de cada tratamiento, entre Hallmarks
    resultados_tratamiento["q_BH"] = bh_adjust(
        resultados_tratamiento["p_empirico"].to_numpy()
    )

    resultados_tratamiento = resultados_tratamiento.sort_values(
        ["q_BH", "p_empirico", "NES"],
        ascending=[True, True, False]
    ).reset_index(drop=True)

    archivo_salida = (
        f"GSEA_HALLMARK_{tratamiento}.csv"
    )

    resultados_tratamiento.to_csv(
        archivo_salida,
        index=False
    )

    n_sig = int(
        (resultados_tratamiento["q_BH"] < 0.05).sum()
    )

    print("\nResultados guardados en:", archivo_salida)
    print("Hallmarks con q_BH < 0.05:", n_sig)

    todos_resultados.append(
        resultados_tratamiento
    )


# ================================================================
# RESULTADOS COMBINADOS
# ================================================================

resultado_final = pd.concat(
    todos_resultados,
    ignore_index=True
)

resultado_final.to_csv(
    "GSEA_HALLMARK_todos_tratamientos.csv",
    index=False
)

significativos = resultado_final[
    resultado_final["q_BH"] < 0.05
].copy()

significativos.to_csv(
    "GSEA_HALLMARK_significativos_q005.csv",
    index=False
)


# ================================================================
# RESUMEN
# ================================================================

resumen = (
    resultado_final
    .assign(
        significativo=lambda x: x["q_BH"] < 0.05,
        positivo=lambda x: (
            (x["q_BH"] < 0.05) &
            (x["NES"] > 0)
        ),
        negativo=lambda x: (
            (x["q_BH"] < 0.05) &
            (x["NES"] < 0)
        )
    )
    .groupby("tratamiento")
    .agg(
        hallmarks_analizados=("hallmark", "size"),
        significativos=("significativo", "sum"),
        positivos=("positivo", "sum"),
        negativos=("negativo", "sum")
    )
    .reset_index()
)

resumen.to_csv(
    "GSEA_HALLMARK_resumen.csv",
    index=False
)


# ================================================================
# METODOLOGÍA
# ================================================================

with open(
    "GSEA_HALLMARK_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "Hallmark preranked enrichment analysis\n"
    )
    f.write(
        "======================================\n\n"
    )
    f.write(
        f"Hallmark file: {GMT_FILE}\n"
    )
    f.write(
        f"Number of permutations per gene set: {N_PERM}\n"
    )
    f.write(
        f"Random seed: {SEED}\n"
    )
    f.write(
        f"Minimum genes per gene set: {MIN_GENES}\n"
    )
    f.write(
        f"Weight exponent: {WEIGHT}\n\n"
    )
    f.write(
        "Genes were ranked by the mean MODZ score across "
        "the four matched primary contexts for each treatment.\n"
    )
    f.write(
        "Enrichment scores were calculated using a weighted "
        "running-sum statistic.\n"
    )
    f.write(
        "Null distributions were generated by randomly sampling "
        "gene sets of identical size from the ranked gene universe.\n"
    )
    f.write(
        "NES values were obtained by normalising the observed ES "
        "against null ES values of the same sign.\n"
    )
    f.write(
        "Empirical one-sided P values were calculated from the "
        "permutation distribution with a +1 correction.\n"
    )
    f.write(
        "Benjamini-Hochberg correction was applied separately "
        "within each treatment across the tested Hallmark sets.\n"
    )


print("\n" + "=" * 80)
print("RESUMEN FINAL")
print("=" * 80)

print(resumen.to_string(index=False))

print("\nArchivos generados:")
print("- GSEA_HALLMARK_estradiol.csv")
print("- GSEA_HALLMARK_calcitriol.csv")
print("- GSEA_HALLMARK_calcipotriol.csv")
print("- GSEA_HALLMARK_todos_tratamientos.csv")
print("- GSEA_HALLMARK_significativos_q005.csv")
print("- GSEA_HALLMARK_resumen.csv")
print("- GSEA_HALLMARK_metodologia.txt")

print("\nFIN DEL GSEA HALLMARK")