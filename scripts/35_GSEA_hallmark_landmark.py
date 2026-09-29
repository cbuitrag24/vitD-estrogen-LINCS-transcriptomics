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

GMT_FILE = "data/h.all.v2026.1.Hs.symbols.gmt"

print("=" * 80)
print("GSEA HALLMARK - ANÁLISIS DE SENSIBILIDAD LANDMARK")
print("=" * 80)

print("\nConfiguración:")
print("Permutaciones por Hallmark:", N_PERM)
print("Semilla:", SEED)
print("Mínimo de genes por Hallmark:", MIN_GENES)
print("Exponente de ponderación:", WEIGHT)


# ================================================================
# FUNCIONES
# ================================================================

def calcular_es(scores, hit_indices, weight=1.0):

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

    pvalues = np.asarray(pvalues, dtype=float)

    n = len(pvalues)
    orden = np.argsort(pvalues)
    p_ordenados = pvalues[orden]

    valores = p_ordenados * n / np.arange(1, n + 1)
    valores = np.minimum.accumulate(valores[::-1])[::-1]
    valores = np.clip(valores, 0, 1)

    q = np.empty(n, dtype=float)
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
        f"No se encontró {GMT_FILE}"
    )

hallmarks = leer_gmt(GMT_FILE)

print("\nHallmarks leídos:", len(hallmarks))

if len(hallmarks) != 50:
    raise ValueError(
        f"Se esperaban 50 Hallmarks y se encontraron {len(hallmarks)}"
    )


# ================================================================
# ANALIZAR CADA TRATAMIENTO
# ================================================================

rng = np.random.default_rng(SEED)

todos_resultados = []

for tratamiento in TRATAMIENTOS:

    print("\n" + "=" * 80)
    print("TRATAMIENTO:", tratamiento.upper())
    print("=" * 80)

    archivo = f"RANKING_{tratamiento}.csv"

    ranking = pd.read_csv(archivo)

    columnas_necesarias = {
        "pr_gene_symbol",
        "mean_modz",
        "pr_is_lm"
    }

    if not columnas_necesarias.issubset(ranking.columns):
        raise ValueError(
            f"{archivo} no contiene las columnas esperadas"
        )

    # ------------------------------------------------------------
    # FILTRAR SOLAMENTE GENES LANDMARK
    # ------------------------------------------------------------

    ranking_landmark = ranking[
        ranking["pr_is_lm"] == 1
    ].copy()

    ranking_landmark["pr_gene_symbol"] = (
        ranking_landmark["pr_gene_symbol"]
        .astype(str)
        .str.strip()
    )

    ranking_landmark["mean_modz"] = pd.to_numeric(
        ranking_landmark["mean_modz"],
        errors="coerce"
    )

    ranking_landmark = ranking_landmark.dropna(
        subset=["pr_gene_symbol", "mean_modz"]
    )

    duplicados = (
        ranking_landmark["pr_gene_symbol"]
        .duplicated()
        .sum()
    )

    if duplicados > 0:
        raise ValueError(
            f"Se encontraron {duplicados} símbolos duplicados "
            f"en los genes landmark de {tratamiento}"
        )

    ranking_landmark = ranking_landmark.sort_values(
        ["mean_modz", "pr_gene_symbol"],
        ascending=[False, True],
        kind="mergesort"
    ).reset_index(drop=True)

    n_total = len(ranking_landmark)

    print("Genes landmark en ranking:", n_total)

    if n_total != 978:
        raise ValueError(
            f"Se esperaban 978 genes landmark, "
            f"pero se encontraron {n_total}"
        )

    genes_ranking = (
        ranking_landmark["pr_gene_symbol"]
        .to_numpy()
    )

    scores = (
        ranking_landmark["mean_modz"]
        .to_numpy(dtype=float)
    )

    mapa_posicion = {
        gen: i
        for i, gen in enumerate(genes_ranking)
    }

    resultados = []

    hallmarks_analizados = 0

    for hallmark, genes_set in hallmarks.items():

        genes_presentes = [
            gen
            for gen in genes_set
            if gen in mapa_posicion
        ]

        n_hits = len(genes_presentes)

        if n_hits < MIN_GENES:
            continue

        hallmarks_analizados += 1

        hit_indices = np.array(
            [mapa_posicion[gen] for gen in genes_presentes],
            dtype=int
        )

        es_obs = calcular_es(
            scores,
            hit_indices,
            weight=WEIGHT
        )

        null_es = np.empty(
            N_PERM,
            dtype=float
        )

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

        if es_obs >= 0:

            null_mismo_signo = null_es[
                null_es >= 0
            ]

            if len(null_mismo_signo) == 0:
                nes = np.nan
            else:
                nes = (
                    es_obs /
                    np.mean(null_mismo_signo)
                )

            p_empirico = (
                np.sum(null_es >= es_obs) + 1
            ) / (N_PERM + 1)

        else:

            null_mismo_signo = null_es[
                null_es < 0
            ]

            if len(null_mismo_signo) == 0:
                nes = np.nan
            else:
                nes = (
                    es_obs /
                    np.mean(
                        np.abs(null_mismo_signo)
                    )
                )

            p_empirico = (
                np.sum(null_es <= es_obs) + 1
            ) / (N_PERM + 1)

        resultados.append({
            "tratamiento": tratamiento,
            "hallmark": hallmark,
            "n_genes_landmark": n_hits,
            "ES": es_obs,
            "NES": nes,
            "p_empirico": p_empirico
        })

    resultados = pd.DataFrame(
        resultados
    )

    if len(resultados) == 0:
        raise ValueError(
            "Ningún Hallmark superó el mínimo "
            "de genes en el espacio landmark"
        )

    resultados["q_BH"] = bh_adjust(
        resultados["p_empirico"].to_numpy()
    )

    resultados = resultados.sort_values(
        ["q_BH", "p_empirico", "NES"],
        ascending=[True, True, False]
    ).reset_index(drop=True)

    resultados.to_csv(
        f"GSEA_LANDMARK_{tratamiento}.csv",
        index=False
    )

    n_sig = int(
        (resultados["q_BH"] < 0.05).sum()
    )

    n_pos = int(
        (
            (resultados["q_BH"] < 0.05) &
            (resultados["NES"] > 0)
        ).sum()
    )

    n_neg = int(
        (
            (resultados["q_BH"] < 0.05) &
            (resultados["NES"] < 0)
        ).sum()
    )

    print("Hallmarks que superan >=15 genes:", hallmarks_analizados)
    print("Hallmarks significativos q<0.05:", n_sig)
    print("  Positivos:", n_pos)
    print("  Negativos:", n_neg)

    todos_resultados.append(
        resultados
    )


# ================================================================
# UNIR LOS TRES TRATAMIENTOS
# ================================================================

resultado_final = pd.concat(
    todos_resultados,
    ignore_index=True
)

resultado_final.to_csv(
    "GSEA_LANDMARK_todos_tratamientos.csv",
    index=False
)

significativos = resultado_final[
    resultado_final["q_BH"] < 0.05
].copy()

significativos.to_csv(
    "GSEA_LANDMARK_significativos_q005.csv",
    index=False
)


# ================================================================
# RESUMEN
# ================================================================

resumen = (
    resultado_final
    .assign(
        significativo=lambda x:
            x["q_BH"] < 0.05,

        positivo=lambda x:
            (x["q_BH"] < 0.05) &
            (x["NES"] > 0),

        negativo=lambda x:
            (x["q_BH"] < 0.05) &
            (x["NES"] < 0)
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
    "GSEA_LANDMARK_resumen.csv",
    index=False
)


# ================================================================
# METODOLOGÍA
# ================================================================

with open(
    "GSEA_LANDMARK_metodologia.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "Hallmark preranked enrichment analysis: "
        "L1000 landmark-gene sensitivity analysis\n"
    )
    f.write(
        "=======================================================\n\n"
    )
    f.write(
        f"Hallmark file: {GMT_FILE}\n"
    )
    f.write(
        "Gene universe: 978 directly measured "
        "L1000 landmark genes\n"
    )
    f.write(
        f"Number of permutations per gene set: {N_PERM}\n"
    )
    f.write(
        f"Random seed: {SEED}\n"
    )
    f.write(
        f"Minimum represented genes per gene set: {MIN_GENES}\n"
    )
    f.write(
        f"Weight exponent: {WEIGHT}\n"
    )
    f.write(
        "Genes were ranked by mean MODZ across "
        "the four matched primary contexts.\n"
    )
    f.write(
        "Empirical one-sided P values were calculated "
        "from random gene sets of identical size with "
        "a +1 correction.\n"
    )
    f.write(
        "Benjamini-Hochberg correction was performed "
        "separately within each treatment.\n"
    )


print("\n" + "=" * 80)
print("RESUMEN FINAL - LANDMARK")
print("=" * 80)

print(
    resumen.to_string(index=False)
)

print("\nArchivos generados:")
print("- GSEA_LANDMARK_estradiol.csv")
print("- GSEA_LANDMARK_calcitriol.csv")
print("- GSEA_LANDMARK_calcipotriol.csv")
print("- GSEA_LANDMARK_todos_tratamientos.csv")
print("- GSEA_LANDMARK_significativos_q005.csv")
print("- GSEA_LANDMARK_resumen.csv")
print("- GSEA_LANDMARK_metodologia.txt")

print("\nFIN DEL GSEA LANDMARK")