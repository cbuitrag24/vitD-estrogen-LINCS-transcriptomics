
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.colors import Normalize
from matplotlib.patches import Patch

# ============================================================
# FIGURE 3 — CORRECTED A4 LANDSCAPE VERSION
# Gene-level convergence between treatments
#
# Size: 29.7 x 21 cm
# PDF vector + PNG 600 dpi
# ============================================================

SOURCE = Path(
    r"C:\Users\petri\LINCS_Estradiol_VitaminaD"
)

REPO = Path(
    r"C:\Users\petri\vitD-estrogen-LINCS-transcriptomics"
)

OUTPUT = REPO / "results" / "figures"

SUMMARY_NAME = "ENRIQUECIMIENTO_resumen_same_direction.csv"
DIRECTION_NAME = "ENRIQUECIMIENTO_solapamiento_resultados.csv"

print("=" * 75)
print("FIGURE 3 — CORRECTED PUBLICATION VERSION")
print("=" * 75)


# ============================================================
# 1. LOCATE ORIGINAL FILES
# ============================================================

def locate_file(filename):
    candidates = []

    for folder in (SOURCE, REPO):
        if folder.is_dir():
            candidates.extend(
                p for p in folder.rglob(filename)
                if p.is_file()
            )

    unique = {
        str(p.resolve()).lower(): p
        for p in candidates
    }

    if len(unique) != 1:
        raise FileNotFoundError(
            f"Expected one copy of {filename}; "
            f"found {len(unique)}:\n"
            + "\n".join(str(p) for p in unique.values())
        )

    return next(iter(unique.values()))


summary_path = locate_file(SUMMARY_NAME)
direction_path = locate_file(DIRECTION_NAME)

print("Summary:", summary_path)
print("Direction:", direction_path)

resumen = pd.read_csv(summary_path)
direccion = pd.read_csv(direction_path)


# ============================================================
# 2. VALIDATE ORIGINAL RESULTS
# ============================================================

required_summary = [
    "analysis_space",
    "comparison",
    "observed_same_direction",
    "expected_same_direction",
    "fold_same_direction"
]

required_direction = [
    "analysis_space",
    "comparison",
    "direction",
    "observed_overlap",
    "expected_overlap",
    "fold_enrichment",
    "p_hypergeom",
    "q_BH"
]

for col in required_summary:
    if col not in resumen.columns:
        raise ValueError(f"Missing summary column: {col}")

for col in required_direction:
    if col not in direccion.columns:
        raise ValueError(f"Missing direction column: {col}")


comparison_order = [
    "calcitriol_vs_calcipotriol",
    "estradiol_vs_calcitriol",
    "estradiol_vs_calcipotriol"
]

spaces = [
    "full_12328_genes",
    "landmark_978_genes"
]

names = {
    "calcitriol_vs_calcipotriol":
        "Calcitriol vs\nCalcipotriol",
    "estradiol_vs_calcitriol":
        "Estradiol vs\nCalcitriol",
    "estradiol_vs_calcipotriol":
        "Estradiol vs\nCalcipotriol"
}

columns = [
    ("full_12328_genes", "up", "Full\nUp"),
    ("full_12328_genes", "down", "Full\nDown"),
    ("landmark_978_genes", "up", "Landmark\nUp"),
    ("landmark_978_genes", "down", "Landmark\nDown")
]

expected_summary = {
    (space, comparison)
    for space in spaces
    for comparison in comparison_order
}

actual_summary = list(zip(
    resumen["analysis_space"],
    resumen["comparison"]
))

if (
    len(actual_summary) != 6
    or len(set(actual_summary)) != 6
    or set(actual_summary) != expected_summary
):
    raise ValueError(
        "Summary must contain exactly six combinations."
    )


expected_direction = {
    (space, comparison, direction)
    for space in spaces
    for comparison in comparison_order
    for direction in ("up", "down")
}

actual_direction = list(zip(
    direccion["analysis_space"],
    direccion["comparison"],
    direccion["direction"]
))

if (
    len(actual_direction) != 12
    or len(set(actual_direction)) != 12
    or set(actual_direction) != expected_direction
):
    raise ValueError(
        "Direction data must contain exactly 12 combinations."
    )


numeric_summary = [
    "observed_same_direction",
    "expected_same_direction",
    "fold_same_direction"
]

numeric_direction = [
    "observed_overlap",
    "expected_overlap",
    "fold_enrichment",
    "p_hypergeom",
    "q_BH"
]

for col in numeric_summary:
    resumen[col] = pd.to_numeric(
        resumen[col], errors="raise"
    )

for col in numeric_direction:
    direccion[col] = pd.to_numeric(
        direccion[col], errors="raise"
    )

if not np.isfinite(
    resumen[numeric_summary].to_numpy()
).all():
    raise ValueError("Invalid summary values.")

if not np.isfinite(
    direccion[numeric_direction].to_numpy()
).all():
    raise ValueError("Invalid direction values.")

if (resumen[numeric_summary] < 0).any().any():
    raise ValueError("Negative summary values.")

if (
    direccion[
        ["observed_overlap", "expected_overlap", "fold_enrichment"]
    ] < 0
).any().any():
    raise ValueError("Negative overlap values.")

if (
    (direccion[["p_hypergeom", "q_BH"]] < 0).any().any()
    or
    (direccion[["p_hypergeom", "q_BH"]] > 1).any().any()
):
    raise ValueError("Invalid p or q values.")


# ============================================================
# 3. PUBLICATION STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white"
})


# ============================================================
# 4. CREATE A4 LANDSCAPE FIGURE
# ============================================================

fig = plt.figure(
    figsize=(29.7 / 2.54, 21 / 2.54),
    facecolor="white"
)

# Fixed axes positions
ax1 = fig.add_axes([0.085, 0.565, 0.385, 0.305])
ax2 = fig.add_axes([0.575, 0.565, 0.385, 0.305])

ax3 = fig.add_axes([0.235, 0.145, 0.635, 0.285])
cax = fig.add_axes([0.89, 0.145, 0.018, 0.285])


# ============================================================
# 5. PANEL LETTERS — SEPARATE FROM TITLES
# ============================================================

fig.text(
    0.025, 0.975, "A",
    fontsize=18,
    fontweight="bold",
    ha="left", va="top"
)

fig.text(
    0.515, 0.975, "B",
    fontsize=18,
    fontweight="bold",
    ha="left", va="top"
)

fig.text(
    0.025, 0.495, "C",
    fontsize=18,
    fontweight="bold",
    ha="left", va="top"
)


# ============================================================
# 6. PANEL A AND B TITLES — FIGURE COORDINATES
# ============================================================

# Titles are NOT attached to the axes.
# This prevents collisions with legends.

fig.text(
    0.2775, 0.945,
    "Full L1000 gene space",
    ha="center",
    va="center",
    fontsize=13,
    fontweight="bold"
)

fig.text(
    0.7675, 0.945,
    "Landmark genes",
    ha="center",
    va="center",
    fontsize=13,
    fontweight="bold"
)


# ============================================================
# 7. OBSERVED / EXPECTED LEGENDS
# ============================================================

observed_color = "#0072B2"
expected_color = "#E69F00"

legend_handles = [
    Patch(
        facecolor=observed_color,
        edgecolor="black",
        linewidth=0.4,
        label="Observed"
    ),
    Patch(
        facecolor=expected_color,
        edgecolor="black",
        linewidth=0.4,
        label="Expected"
    )
]

# Legend A — below its title
fig.legend(
    handles=legend_handles,
    loc="center",
    bbox_to_anchor=(0.2775, 0.905),
    ncol=2,
    frameon=False,
    fontsize=11,
    handlelength=1.5,
    columnspacing=2
)

# Legend B — below its title
fig.legend(
    handles=legend_handles,
    loc="center",
    bbox_to_anchor=(0.7675, 0.905),
    ncol=2,
    frameon=False,
    fontsize=11,
    handlelength=1.5,
    columnspacing=2
)


# ============================================================
# 8. PLOT OBSERVED VS EXPECTED
# ============================================================

def plot_observed_expected(ax, space):

    data = (
        resumen.loc[
            resumen["analysis_space"] == space
        ]
        .set_index("comparison")
        .loc[comparison_order]
    )

    x = np.arange(3)
    width = 0.32

    observed = data[
        "observed_same_direction"
    ].to_numpy()

    expected = data[
        "expected_same_direction"
    ].to_numpy()

    fold = data[
        "fold_same_direction"
    ].to_numpy()

    bars_observed = ax.bar(
        x - width / 2,
        observed,
        width,
        color=observed_color,
        edgecolor="black",
        linewidth=0.4,
        zorder=3
    )

    bars_expected = ax.bar(
        x + width / 2,
        expected,
        width,
        color=expected_color,
        edgecolor="black",
        linewidth=0.4,
        zorder=3
    )

    max_height = float(
        max(np.max(observed), np.max(expected))
    )

    ax.set_ylim(0, max_height * 1.34)
    ax.set_xlim(-0.52, 2.52)

    ax.set_xticks(x)

    ax.set_xticklabels(
        [names[c] for c in comparison_order],
        fontsize=9.4
    )

    ax.set_ylabel("Same-direction genes")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        direction="out",
        length=4,
        width=0.8
    )

    for bar, value in zip(bars_observed, observed):
        ax.annotate(
            f"{value:.0f}",
            (
                bar.get_x() + bar.get_width() / 2,
                bar.get_height()
            ),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            color="#202020"
        )

    for bar, value in zip(bars_expected, expected):
        ax.annotate(
            f"{value:.1f}",
            (
                bar.get_x() + bar.get_width() / 2,
                bar.get_height()
            ),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            color="#202020"
        )

    for i, value in enumerate(fold):
        height = max(observed[i], expected[i])

        ax.text(
            i,
            height + max_height * 0.14,
            f"{value:.2f}×",
            ha="center",
            va="center",
            fontsize=11.5,
            fontweight="bold",
            color="#202020"
        )

    ax.grid(
        axis="y",
        color="#E6E6E6",
        linewidth=0.6,
        zorder=0
    )


plot_observed_expected(
    ax1,
    "full_12328_genes"
)

plot_observed_expected(
    ax2,
    "landmark_978_genes"
)


# ============================================================
# 9. PANEL C — DIRECTION-SPECIFIC ENRICHMENT
# ============================================================

fold_matrix = np.empty((3, 4))
q_matrix = np.empty((3, 4))

for i, comparison in enumerate(comparison_order):

    for j, (space, direction, _) in enumerate(columns):

        row = direccion.loc[
            (direccion["comparison"] == comparison)
            &
            (direccion["analysis_space"] == space)
            &
            (direccion["direction"] == direction)
        ]

        if len(row) != 1:
            raise ValueError(
                f"Missing or duplicated combination: "
                f"{comparison}, {space}, {direction}"
            )

        fold_matrix[i, j] = (
            row["fold_enrichment"].iloc[0]
        )

        q_matrix[i, j] = (
            row["q_BH"].iloc[0]
        )


norm = Normalize(
    vmin=float(np.min(fold_matrix)),
    vmax=float(np.max(fold_matrix))
)

if norm.vmin == norm.vmax:
    raise ValueError(
        "All enrichment values are identical."
    )

im = ax3.imshow(
    fold_matrix,
    cmap="viridis",
    norm=norm,
    aspect="auto",
    interpolation="nearest"
)

ax3.set_xticks(range(4))
ax3.set_xticklabels(
    [col[2] for col in columns],
    fontsize=11
)

ax3.set_yticks(range(3))
ax3.set_yticklabels(
    [
        names[c].replace("\n", " ")
        for c in comparison_order
    ],
    fontsize=10.5
)

ax3.set_title(
    "Direction-specific overlap enrichment",
    fontsize=13,
    fontweight="bold",
    pad=13
)

ax3.tick_params(length=0, pad=8)

# White borders between cells
ax3.set_xticks(
    np.arange(-0.5, 4, 1),
    minor=True
)

ax3.set_yticks(
    np.arange(-0.5, 3, 1),
    minor=True
)

ax3.grid(
    which="minor",
    color="white",
    linewidth=2
)

ax3.tick_params(
    which="minor",
    bottom=False,
    left=False
)


# ============================================================
# 10. HEATMAP LABELS AND SIGNIFICANCE
# ============================================================

cmap = plt.get_cmap("viridis")

for i in range(3):
    for j in range(4):

        fold = fold_matrix[i, j]
        q = q_matrix[i, j]

        if q < 0.001:
            stars = "***"
        elif q < 0.01:
            stars = "**"
        elif q < 0.05:
            stars = "*"
        else:
            stars = ""

        rgba = cmap(norm(fold))

        luminance = (
            0.2126 * rgba[0]
            + 0.7152 * rgba[1]
            + 0.0722 * rgba[2]
        )

        text_color = (
            "#111111"
            if luminance > 0.52
            else "white"
        )

        ax3.text(
            j,
            i,
            f"{fold:.2f}×{stars}",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold" if stars else "normal",
            color=text_color
        )


# ============================================================
# 11. COLOUR BAR AND STATISTICAL NOTE
# ============================================================

cbar = fig.colorbar(im, cax=cax)

cbar.set_label(
    "Fold enrichment",
    fontsize=11,
    labelpad=12
)

cbar.ax.tick_params(labelsize=10)

fig.text(
    0.55,
    0.065,
    "* BH-adjusted P < 0.05; "
    "** P < 0.01; "
    "*** P < 0.001",
    ha="center",
    va="center",
    fontsize=10,
    color="#333333"
)


# ============================================================
# 12. EXPORT — OVERWRITE PREVIOUS FILES
# ============================================================

OUTPUT.mkdir(parents=True, exist_ok=True)

pdf_path = (
    OUTPUT / "Figure_3_gene_level_convergence.pdf"
)

png_path = (
    OUTPUT / "Figure_3_gene_level_convergence.png"
)

csv_path = (
    OUTPUT / "Figure_3_gene_overlap_values.csv"
)

# Exact A4 dimensions: no tight bounding box
fig.savefig(
    pdf_path,
    format="pdf",
    bbox_inches=None,
    facecolor="white"
)

fig.savefig(
    png_path,
    format="png",
    dpi=600,
    bbox_inches=None,
    facecolor="white"
)

audit_table = direccion[
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

audit_table.to_csv(
    csv_path,
    index=False
)

plt.close(fig)

print("\n" + "=" * 75)
print("FIGURE 3 GENERATED SUCCESSFULLY")
print("=" * 75)
print("Size: A4 landscape — 29.7 x 21 cm")
print("PDF:", pdf_path)
print("PNG:", png_path)
print("Audit CSV:", csv_path)
