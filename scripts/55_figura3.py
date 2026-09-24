import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

print("=" * 90)
print("FIGURA 3 — CONVERGENCIA GÉNICA ENTRE TRATAMIENTOS")
print("=" * 90)

# ============================================================
# ARCHIVOS
# ============================================================

archivo_resumen = "ENRIQUECIMIENTO_resumen_same_direction.csv"
archivo_direccion = "ENRIQUECIMIENTO_solapamiento_resultados.csv"

resumen = pd.read_csv(archivo_resumen)
direccion = pd.read_csv(archivo_direccion)

print("\nResumen same-direction:")
print(resumen.to_string(index=False))

print("\nResultados por dirección:")
print(direccion.to_string(index=False))

# ============================================================
# NOMBRES PARA LA FIGURA
# ============================================================

mapa_nombres = {
    "calcitriol_vs_calcipotriol":
        "Calcitriol vs\nCalcipotriol",

    "estradiol_vs_calcitriol":
        "Estradiol vs\nCalcitriol",

    "estradiol_vs_calcipotriol":
        "Estradiol vs\nCalcipotriol"
}

orden = [
    "calcitriol_vs_calcipotriol",
    "estradiol_vs_calcitriol",
    "estradiol_vs_calcipotriol"
]

# ============================================================
# FIGURA
# ============================================================

fig = plt.figure(
    figsize=(12.5, 9.0),
    constrained_layout=True
)

gs = fig.add_gridspec(
    2,
    2,
    height_ratios=[1, 1.08]
)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])

# ============================================================
# FUNCIÓN PARA PANELES A Y B
# ============================================================

def plot_observed_expected(ax, data, title, panel_letter):

    data = (
        data
        .set_index("comparison")
        .loc[orden]
        .reset_index()
    )

    x = np.arange(len(data))
    width = 0.34

    observed = data["observed_same_direction"].values
    expected = data["expected_same_direction"].values
    fold = data["fold_same_direction"].values

    bars_obs = ax.bar(
        x - width / 2,
        observed,
        width,
        label="Observed"
    )

    bars_exp = ax.bar(
        x + width / 2,
        expected,
        width,
        label="Expected"
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        [
            mapa_nombres[c]
            for c in data["comparison"]
        ],
        fontsize=9
    )

    ax.set_ylabel(
        "Same-direction genes"
    )

    ax.set_title(
        title,
        fontsize=11
    )

    ax.legend(
        frameon=False,
        fontsize=9
    )

    ax.text(
        -0.11,
        1.05,
        panel_letter,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top"
    )

    # Etiquetas observadas
    for bar, value in zip(
        bars_obs,
        observed
    ):

        ax.text(
            bar.get_x()
            + bar.get_width() / 2,
            bar.get_height()
            + max(observed) * 0.025,
            f"{int(value)}",
            ha="center",
            va="bottom",
            fontsize=8.5
        )

    # Etiquetas esperadas
    for bar, value in zip(
        bars_exp,
        expected
    ):

        ax.text(
            bar.get_x()
            + bar.get_width() / 2,
            bar.get_height()
            + max(observed) * 0.025,
            f"{value:.1f}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    # Fold enrichment
    for i, f in enumerate(fold):

        ymax = max(
            observed[i],
            expected[i]
        )

        ax.text(
            i,
            ymax
            + max(observed) * 0.11,
            f"{f:.2f}×",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold"
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ============================================================
# PANEL A — FULL
# ============================================================

full = resumen[
    resumen["analysis_space"]
    == "full_12328_genes"
].copy()

plot_observed_expected(
    ax1,
    full,
    "Full L1000 gene space",
    "A"
)

# ============================================================
# PANEL B — LANDMARK
# ============================================================

landmark = resumen[
    resumen["analysis_space"]
    == "landmark_978_genes"
].copy()

plot_observed_expected(
    ax2,
    landmark,
    "Landmark genes",
    "B"
)

# ============================================================
# PANEL C — FOLD ENRICHMENT POR DIRECCIÓN
# ============================================================

ax3.text(
    -0.055,
    1.05,
    "C",
    transform=ax3.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

# Orden de columnas
columnas = [
    ("full_12328_genes", "up", "Full\nUp"),
    ("full_12328_genes", "down", "Full\nDown"),
    ("landmark_978_genes", "up", "Landmark\nUp"),
    ("landmark_978_genes", "down", "Landmark\nDown")
]

matriz_fold = np.zeros(
    (
        len(orden),
        len(columnas)
    )
)

matriz_q = np.ones(
    (
        len(orden),
        len(columnas)
    )
)

for i, comparacion in enumerate(orden):

    for j, (
        espacio,
        direc,
        _
    ) in enumerate(columnas):

        fila = direccion[
            (
                direccion["comparison"]
                == comparacion
            )
            &
            (
                direccion["analysis_space"]
                == espacio
            )
            &
            (
                direccion["direction"]
                == direc
            )
        ]

        if len(fila) != 1:
            raise ValueError(
                f"Problema al localizar: "
                f"{comparacion}, {espacio}, {direc}"
            )

        matriz_fold[i, j] = (
            fila["fold_enrichment"]
            .iloc[0]
        )

        matriz_q[i, j] = (
            fila["q_BH"]
            .iloc[0]
        )

# Heatmap
im = ax3.imshow(
    matriz_fold,
    aspect="auto",
    cmap="viridis"
)

ax3.set_xticks(
    np.arange(
        len(columnas)
    )
)

ax3.set_xticklabels(
    [
        c[2]
        for c in columnas
    ],
    fontsize=9
)

ax3.set_yticks(
    np.arange(
        len(orden)
    )
)

ax3.set_yticklabels(
    [
        mapa_nombres[c]
        .replace("\n", " ")
        for c in orden
    ],
    fontsize=9
)

ax3.set_title(
    "Direction-specific overlap enrichment",
    fontsize=11
)

# Texto dentro de cada celda
for i in range(
    matriz_fold.shape[0]
):

    for j in range(
        matriz_fold.shape[1]
    ):

        fold = matriz_fold[i, j]
        q = matriz_q[i, j]

        if q < 0.001:
            estrella = "***"
        elif q < 0.01:
            estrella = "**"
        elif q < 0.05:
            estrella = "*"
        else:
            estrella = ""

        ax3.text(
            j,
            i,
            f"{fold:.2f}×{estrella}",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold"
            if q < 0.05
            else "normal",
            color="white"
            if fold > 1.7
            else "black"
        )

# Colorbar
cbar = fig.colorbar(
    im,
    ax=ax3,
    fraction=0.025,
    pad=0.02
)

cbar.set_label(
    "Fold enrichment",
    fontsize=9
)

# Nota estadística
ax3.text(
    0.5,
    -0.22,
    "* BH-adjusted P < 0.05; ** P < 0.01; *** P < 0.001",
    transform=ax3.transAxes,
    ha="center",
    va="top",
    fontsize=8.5
)

# ============================================================
# GUARDAR
# ============================================================

fig.savefig(
    "Figure_3_gene_level_convergence.pdf",
    bbox_inches="tight"
)

fig.savefig(
    "Figure_3_gene_level_convergence.png",
    dpi=600,
    bbox_inches="tight"
)

# Tabla auditiva utilizada
tabla_figura3 = direccion[
    [
        "analysis_space",
        "comparison",
        "direction",
        "observed_overlap",
        "expected_overlap",
        "fold_enrichment",
        "p_hypergeom",
        "q_BH"
    ]
].copy()

tabla_figura3.to_csv(
    "Figure_3_gene_overlap_values.csv",
    index=False
)

plt.close(fig)

print("\n" + "=" * 90)
print("ARCHIVOS GENERADOS")
print("=" * 90)

print(
    "- Figure_3_gene_level_convergence.pdf"
)

print(
    "- Figure_3_gene_level_convergence.png"
)

print(
    "- Figure_3_gene_overlap_values.csv"
)

print("\nFIGURA 3 COMPLETADA")