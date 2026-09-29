import pandas as pd
import numpy as np
from pathlib import Path

print("=" * 90)
print("LEADING-EDGE GENES DEL NÚCLEO HALLMARK ROBUSTO")
print("=" * 90)

GMT_FILE = "data/h.all.v2026.1.Hs.symbols.gmt"

TRATAMIENTOS = {
    "estradiol": "RANKING_estradiol.csv",
    "calcitriol": "RANKING_calcitriol.csv",
    "calcipotriol": "RANKING_calcipotriol.csv",
}

HALLMARKS_NUCLEO = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
]

# ------------------------------------------------------------
# 1. LEER GMT
# ------------------------------------------------------------

def leer_gmt(archivo):
    conjuntos = {}

    with open(archivo, "r", encoding="utf-8") as f:
        for linea in f:
            partes = linea.rstrip("\n").split("\t")

            if len(partes) >= 3:
                nombre = partes[0]
                genes = set(partes[2:])
                conjuntos[nombre] = genes

    return conjuntos


hallmarks = leer_gmt(GMT_FILE)

print(f"\nHallmarks leídos: {len(hallmarks)}")

for h in HALLMARKS_NUCLEO:
    if h not in hallmarks:
        raise ValueError(f"No se encontró {h} en el GMT")

# ------------------------------------------------------------
# 2. FUNCIÓN PARA CALCULAR ES Y LEADING EDGE
# ------------------------------------------------------------

def calcular_leading_edge(ranking, genes_set):

    ranking = ranking.copy()

    ranking = ranking.dropna(
        subset=["pr_gene_symbol", "mean_modz"]
    )

    ranking = (
        ranking
        .sort_values("mean_modz", ascending=False)
        .drop_duplicates("pr_gene_symbol")
        .reset_index(drop=True)
    )

    genes_ranking = ranking["pr_gene_symbol"].astype(str).values
    scores = ranking["mean_modz"].astype(float).values

    hits = np.array(
        [g in genes_set for g in genes_ranking],
        dtype=bool
    )

    Nh = hits.sum()
    N = len(ranking)

    if Nh == 0:
        return None

    if Nh == N:
        raise ValueError(
            "El Hallmark ocupa todo el ranking; ES no definido."
        )

    # Peso GSEA p=1
    hit_weights = np.abs(scores[hits])

    if hit_weights.sum() == 0:
        hit_weights = np.ones_like(hit_weights)

    hit_weights = hit_weights / hit_weights.sum()

    running = np.zeros(N, dtype=float)

    hit_counter = 0
    miss_penalty = 1.0 / (N - Nh)

    acumulado = 0.0

    for i in range(N):

        if hits[i]:
            acumulado += hit_weights[hit_counter]
            hit_counter += 1
        else:
            acumulado -= miss_penalty

        running[i] = acumulado

    max_es = running.max()
    min_es = running.min()

    if abs(max_es) >= abs(min_es):
        ES = max_es
        pico = int(np.argmax(running))
        direccion = "positivo"

        # Para enriquecimiento positivo:
        # genes del set encontrados hasta el máximo ES.
        mascara_le = hits.copy()
        mascara_le[(pico + 1):] = False

    else:
        ES = min_es
        pico = int(np.argmin(running))
        direccion = "negativo"

        # Para enriquecimiento negativo:
        # genes del set desde el mínimo ES hacia el final.
        mascara_le = hits.copy()
        mascara_le[:pico] = False

    indices_le = np.where(mascara_le)[0]

    leading = ranking.iloc[indices_le].copy()

    leading["rank"] = indices_le + 1
    leading["ES"] = ES
    leading["direccion_ES"] = direccion

    return {
        "ES": ES,
        "direccion": direccion,
        "pico_rank": pico + 1,
        "n_set_en_ranking": int(Nh),
        "n_leading_edge": len(leading),
        "leading": leading
    }


# ------------------------------------------------------------
# 3. ANALIZAR FULL Y LANDMARK
# ------------------------------------------------------------

todos_leading = []
resumen = []

for tratamiento, archivo in TRATAMIENTOS.items():

    print("\n" + "=" * 90)
    print(f"TRATAMIENTO: {tratamiento.upper()}")
    print("=" * 90)

    ranking_total = pd.read_csv(archivo)

    columnas_necesarias = {
        "pr_gene_symbol",
        "mean_modz",
        "pr_is_lm"
    }

    if not columnas_necesarias.issubset(ranking_total.columns):
        faltantes = columnas_necesarias - set(ranking_total.columns)
        raise ValueError(
            f"Faltan columnas en {archivo}: {faltantes}"
        )

    espacios = {
        "full": ranking_total.copy(),
        "landmark": ranking_total[
            ranking_total["pr_is_lm"] == 1
        ].copy()
    }

    for espacio, ranking in espacios.items():

        print(f"\nEspacio: {espacio.upper()}")
        print(f"Genes en ranking: {len(ranking)}")

        for hallmark in HALLMARKS_NUCLEO:

            resultado = calcular_leading_edge(
                ranking,
                hallmarks[hallmark]
            )

            if resultado is None:
                print(f"{hallmark}: sin genes")
                continue

            print(
                f"{hallmark}: "
                f"ES={resultado['ES']:.4f} | "
                f"dirección={resultado['direccion']} | "
                f"genes_set={resultado['n_set_en_ranking']} | "
                f"leading_edge={resultado['n_leading_edge']}"
            )

            leading = resultado["leading"].copy()

            leading["tratamiento"] = tratamiento
            leading["espacio"] = espacio
            leading["hallmark"] = hallmark

            todos_leading.append(leading)

            resumen.append({
                "tratamiento": tratamiento,
                "espacio": espacio,
                "hallmark": hallmark,
                "ES": resultado["ES"],
                "direccion": resultado["direccion"],
                "pico_rank": resultado["pico_rank"],
                "genes_set_en_ranking":
                    resultado["n_set_en_ranking"],
                "n_leading_edge":
                    resultado["n_leading_edge"]
            })

# ------------------------------------------------------------
# 4. UNIR Y GUARDAR
# ------------------------------------------------------------

leading_df = pd.concat(
    todos_leading,
    ignore_index=True
)

resumen_df = pd.DataFrame(resumen)

leading_df.to_csv(
    "LEADING_EDGE_nucleo_robusto_todos.csv",
    index=False
)

resumen_df.to_csv(
    "LEADING_EDGE_nucleo_robusto_resumen.csv",
    index=False
)

# ------------------------------------------------------------
# 5. GENES COMPARTIDOS ENTRE LOS TRES TRATAMIENTOS
# ------------------------------------------------------------

compartidos_resumen = []
compartidos_genes = []

for espacio in ["full", "landmark"]:

    print("\n" + "=" * 90)
    print(f"GENES LEADING-EDGE COMPARTIDOS - {espacio.upper()}")
    print("=" * 90)

    for hallmark in HALLMARKS_NUCLEO:

        genes_por_tratamiento = {}

        for tratamiento in TRATAMIENTOS:

            sub = leading_df[
                (leading_df["espacio"] == espacio)
                &
                (leading_df["hallmark"] == hallmark)
                &
                (leading_df["tratamiento"] == tratamiento)
            ]

            genes_por_tratamiento[tratamiento] = set(
                sub["pr_gene_symbol"].astype(str)
            )

        comunes_3 = (
            genes_por_tratamiento["estradiol"]
            &
            genes_por_tratamiento["calcitriol"]
            &
            genes_por_tratamiento["calcipotriol"]
        )

        e2_calcitriol = (
            genes_por_tratamiento["estradiol"]
            &
            genes_por_tratamiento["calcitriol"]
        )

        e2_calcipotriol = (
            genes_por_tratamiento["estradiol"]
            &
            genes_por_tratamiento["calcipotriol"]
        )

        calcitriol_calcipotriol = (
            genes_por_tratamiento["calcitriol"]
            &
            genes_por_tratamiento["calcipotriol"]
        )

        print(f"\n{hallmark}")
        print(
            f"  Compartidos por los 3 tratamientos: "
            f"{len(comunes_3)}"
        )
        print(
            f"  Estradiol–calcitriol: "
            f"{len(e2_calcitriol)}"
        )
        print(
            f"  Estradiol–calcipotriol: "
            f"{len(e2_calcipotriol)}"
        )
        print(
            f"  Calcitriol–calcipotriol: "
            f"{len(calcitriol_calcipotriol)}"
        )

        if comunes_3:
            print(
                "  Genes comunes 3: "
                + ", ".join(sorted(comunes_3))
            )

        compartidos_resumen.append({
            "espacio": espacio,
            "hallmark": hallmark,
            "n_comunes_3": len(comunes_3),
            "n_estradiol_calcitriol":
                len(e2_calcitriol),
            "n_estradiol_calcipotriol":
                len(e2_calcipotriol),
            "n_calcitriol_calcipotriol":
                len(calcitriol_calcipotriol)
        })

        for gen in sorted(comunes_3):
            compartidos_genes.append({
                "espacio": espacio,
                "hallmark": hallmark,
                "gene": gen
            })

# ------------------------------------------------------------
# 6. GUARDAR GENES COMPARTIDOS
# ------------------------------------------------------------

pd.DataFrame(compartidos_resumen).to_csv(
    "LEADING_EDGE_compartidos_resumen.csv",
    index=False
)

pd.DataFrame(compartidos_genes).to_csv(
    "LEADING_EDGE_genes_comunes_3.csv",
    index=False
)

# ------------------------------------------------------------
# 7. GENES REPETIDOS ENTRE LAS 3 VÍAS ROBUSTAS
# ------------------------------------------------------------

conteo = (
    pd.DataFrame(compartidos_genes)
    .groupby(["espacio", "gene"])
    .size()
    .reset_index(name="n_hallmarks")
)

conteo = conteo[
    conteo["n_hallmarks"] >= 2
].sort_values(
    ["espacio", "n_hallmarks", "gene"],
    ascending=[True, False, True]
)

conteo.to_csv(
    "LEADING_EDGE_genes_recurrentes_nucleo.csv",
    index=False
)

print("\n" + "=" * 90)
print("GENES COMPARTIDOS POR LOS 3 TRATAMIENTOS")
print("Y PRESENTES EN >=2 HALLMARKS DEL NÚCLEO")
print("=" * 90)

if len(conteo) == 0:
    print("No se encontraron genes recurrentes.")
else:
    print(conteo.to_string(index=False))

print("\nArchivos generados:")
print("- LEADING_EDGE_nucleo_robusto_todos.csv")
print("- LEADING_EDGE_nucleo_robusto_resumen.csv")
print("- LEADING_EDGE_compartidos_resumen.csv")
print("- LEADING_EDGE_genes_comunes_3.csv")
print("- LEADING_EDGE_genes_recurrentes_nucleo.csv")

print("\nFIN")