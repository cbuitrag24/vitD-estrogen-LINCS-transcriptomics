import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

print("=" * 90)
print("FIGURA 2 — VARIACIÓN GLOBAL Y SIMILITUD ENTRE TRATAMIENTOS")
print("=" * 90)

# ============================================================
# PANEL A — PERMANOVA
# ============================================================

permanova = pd.DataFrame({
    "Gene space": [
        "Full",
        "Full",
        "Landmark",
        "Landmark"
    ],
    "Source": [
        "Treatment",
        "Context / block",
        "Treatment",
        "Context / block"
    ],
    "R2": [
        0.1380151871065361,
        0.4252446816912371,
        0.14877538546041968,
        0.42440058382449297
    ]
})

print("\nPERMANOVA:")
print(permanova.to_string(index=False))

# p exactos del efecto de tratamiento
p_full = 0.5740740740740741
p_landmark = 0.2916666666666667

# ============================================================
# PANEL B — SIMILITUD POR CONTEXTO
# ============================================================

archivo_sim = Path("ESTADISTICA_similitud_por_contexto.csv")

if not archivo_sim.exists():
    raise FileNotFoundError(
        "No se encontró ESTADISTICA_similitud_por_contexto.csv"
    )

sim = pd.read_csv(archivo_sim)

print("\nArchivo de similitud encontrado.")
print("Filas:", len(sim))
print("Columnas:")
print(list(sim.columns))

# ------------------------------------------------------------
# Verificar estructura real del CSV
# ------------------------------------------------------------

columnas_necesarias = [
    "analysis",
    "metric",
    "cell_id",
    "vitd_vitd",
    "e2_calcitriol",
    "e2_calcipotriol"
]

faltantes = [
    c for c in columnas_necesarias
    if c not in sim.columns
]

if faltantes:
    raise ValueError(
        f"Faltan columnas necesarias en el CSV: {faltantes}"
    )

# ------------------------------------------------------------
# Usar sólo Full + Pearson
# ------------------------------------------------------------

sim_filtrado = sim[
    (
        sim["analysis"]
        .astype(str)
        .str.lower()
        .str.contains("full")
    )
    &
    (
        sim["metric"]
        .astype(str)
        .str.lower()
        .str.contains("pearson")
    )
].copy()

print("\nFilas Full + Pearson:", len(sim_filtrado))

if len(sim_filtrado) == 0:
    raise ValueError(
        "No se encontraron filas correspondientes a Full + Pearson."
    )

# ------------------------------------------------------------
# Convertir formato ancho a largo
# ------------------------------------------------------------

sim_plot = sim_filtrado.melt(
    id_vars=["cell_id"],
    value_vars=[
        "vitd_vitd",
        "e2_calcitriol",
        "e2_calcipotriol"
    ],
    var_name="comparison_code",
    value_name="pearson"
)

mapa_comparaciones = {
    "vitd_vitd":
        "Calcitriol vs Calcipotriol",

    "e2_calcitriol":
        "Estradiol vs Calcitriol",

    "e2_calcipotriol":
        "Estradiol vs Calcipotriol"
}

sim_plot["comparison_label"] = (
    sim_plot["comparison_code"]
    .map(mapa_comparaciones)
)

orden_comparaciones = [
    "Calcitriol vs Calcipotriol",
    "Estradiol vs Calcitriol",
    "Estradiol vs Calcipotriol"
]

print("\nDatos utilizados para Panel B:")
print(
    sim_plot[
        [
            "cell_id",
            "comparison_label",
            "pearson"
        ]
    ].to_string(
        index=False
    )
)

# Verificar que estén las tres comparaciones
presentes = set(
    sim_plot["comparison_label"]
)

faltantes = [
    x for x in orden_comparaciones
    if x not in presentes
]

if faltantes:
    raise ValueError(
        f"Faltan comparaciones: {faltantes}"
    )

# ============================================================
# CREAR FIGURA
# ============================================================

fig = plt.figure(
    figsize=(11.5, 5.3),
    constrained_layout=True
)

gs = fig.add_gridspec(
    1,
    2,
    width_ratios=[0.9, 1.35]
)

# ============================================================
# PANEL A — PERMANOVA R²
# ============================================================

ax1 = fig.add_subplot(gs[0, 0])

ax1.text(
    -0.14,
    1.04,
    "A",
    transform=ax1.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

espacios = [
    "Full",
    "Landmark"
]

x = np.arange(len(espacios))

width = 0.32

r2_treatment = [
    0.1380151871065361,
    0.14877538546041968
]

r2_context = [
    0.4252446816912371,
    0.42440058382449297
]

bars1 = ax1.bar(
    x - width / 2,
    r2_treatment,
    width,
    label="Treatment"
)

bars2 = ax1.bar(
    x + width / 2,
    r2_context,
    width,
    label="Context / block"
)

ax1.set_xticks(x)
ax1.set_xticklabels(espacios)

ax1.set_ylabel(
    "Variance explained (R²)"
)

ax1.set_ylim(
    0,
    0.52
)

ax1.set_title(
    "Global variance partitioning",
    fontsize=11
)

ax1.legend(
    frameon=False,
    fontsize=9
)

for bars in [bars1, bars2]:

    for bar in bars:

        valor = bar.get_height()

        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            valor + 0.012,
            f"{valor:.2f}",
            ha="center",
            va="bottom",
            fontsize=8.5
        )

# p exactos sólo para Treatment
ax1.text(
    x[0] - width / 2,
    r2_treatment[0] + 0.065,
    f"p = {p_full:.3f}",
    ha="center",
    fontsize=8
)

ax1.text(
    x[1] - width / 2,
    r2_treatment[1] + 0.065,
    f"p = {p_landmark:.3f}",
    ha="center",
    fontsize=8
)

ax1.text(
    0.5,
    -0.17,
    "Exact within-block permutations for treatment;\n"
    "context/block R² shown descriptively.",
    transform=ax1.transAxes,
    ha="center",
    va="top",
    fontsize=8
)

# ============================================================
# PANEL B — PEARSON POR CONTEXTO
# ============================================================

ax2 = fig.add_subplot(gs[0, 1])

ax2.text(
    -0.10,
    1.04,
    "B",
    transform=ax2.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

contextos = [
    "A549",
    "HA1E",
    "PC3",
    "VCAP"
]

offsets = {
    "A549": -0.12,
    "HA1E": -0.04,
    "PC3": 0.04,
    "VCAP": 0.12
}

markers = {
    "A549": "o",
    "HA1E": "s",
    "PC3": "^",
    "VCAP": "D"
}

for i, comparacion in enumerate(
    orden_comparaciones
):

    temp = sim_plot[
        sim_plot["comparison_label"]
        == comparacion
    ].copy()

    for contexto in contextos:

        fila = temp[
            temp["cell_id"].astype(str)
            == contexto
        ]

        if len(fila) == 0:
            continue

        valor = float(
            fila["pearson"].iloc[0]
        )

        ax2.scatter(
            i + offsets[contexto],
            valor,
            s=58,
            marker=markers[contexto],
            label=contexto
            if i == 0
            else None,
            edgecolors="black",
            linewidths=0.5,
            zorder=3
        )

    valores = pd.to_numeric(
        temp["pearson"],
        errors="coerce"
    ).dropna()

    media = valores.mean()

    ax2.plot(
        [
            i - 0.23,
            i + 0.23
        ],
        [
            media,
            media
        ],
        linewidth=2.2, 
        color="black",
        zorder=2
    )

    ax2.text(
        i,
        media + 0.035,
        f"mean = {media:.2f}",
        ha="center",
        va="bottom",
        fontsize=8
    )

ax2.axhline(
    0,
    linewidth=0.6
)

ax2.set_xticks(
    np.arange(
        len(orden_comparaciones)
    )
)

ax2.set_xticklabels(
    [
        "Calcitriol vs\nCalcipotriol",
        "Estradiol vs\nCalcitriol",
        "Estradiol vs\nCalcipotriol"
    ],
    fontsize=8.5
)

ax2.set_ylabel(
    "Pearson correlation"
)

ax2.set_title(
    "Pairwise transcriptional similarity across contexts",
    fontsize=11
)

ax2.legend(
    title="Cell line",
    frameon=False,
    fontsize=8,
    title_fontsize=8,
    loc="best"
)

# ============================================================
# GUARDAR
# ============================================================

fig.savefig(
    "Figure_2_global_variance_similarity.pdf",
    bbox_inches="tight"
)

fig.savefig(
    "Figure_2_global_variance_similarity.png",
    dpi=600,
    bbox_inches="tight"
)

permanova.to_csv(
    "Figure_2_PERMANOVA_values.csv",
    index=False
)

plt.close(fig)

print("\n" + "=" * 90)
print("ARCHIVOS GENERADOS")
print("=" * 90)

print(
    "- Figure_2_global_variance_similarity.pdf"
)

print(
    "- Figure_2_global_variance_similarity.png"
)

print(
    "- Figure_2_PERMANOVA_values.csv"
)

print("\nFIGURA 2 COMPLETADA")