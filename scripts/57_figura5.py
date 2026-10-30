
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.colors import TwoSlopeNorm

# ============================================================
# FIGURE 5 — CELLULAR CONTEXT DEPENDENCE
# Publication version
# A4 landscape: 29.7 x 21 cm
# Vector PDF + PNG 600 dpi
# ============================================================

SOURCE = Path(
    r"C:\Users\petri\LINCS_Estradiol_VitaminaD"
)

REPO = Path(
    r"C:\Users\petri\vitD-estrogen-LINCS-transcriptomics"
)

OUTPUT = REPO / "results" / "figures"

print("=" * 75)
print("FIGURE 5 — CELLULAR CONTEXT DEPENDENCE")
print("=" * 75)


# ============================================================
# 1. LOCATE ORIGINAL DATA FILES
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


files = {
    "full": "HALLMARK_CONTEXTO_matriz.csv",
    "landmark": "HALLMARK_CONTEXTO_LM_consenso.csv",
    "summary_full": "HALLMARK_CONTEXTO_resumen.csv",
    "summary_landmark": "HALLMARK_CONTEXTO_LM_resumen.csv"
}

paths = {
    key: locate_file(filename)
    for key, filename in files.items()
}

for key, path in paths.items():
    print(f"{key}: {path}")

full = pd.read_csv(paths["full"])
landmark = pd.read_csv(paths["landmark"])

summary_full = pd.read_csv(paths["summary_full"])
summary_landmark = pd.read_csv(paths["summary_landmark"])


# ============================================================
# 2. DEFINE ORIGINAL ORDER
# ============================================================

contexts = [
    "A549",
    "HA1E",
    "PC3",
    "VCAP"
]

hallmarks = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1"
]

treatments = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

hallmark_labels = {
    "HALLMARK_E2F_TARGETS": "E2F TARGETS",
    "HALLMARK_G2M_CHECKPOINT": "G2M CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1": "MYC TARGETS V1"
}

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol"
}

landmark_columns = {
    "estradiol": "estradiol_LM",
    "calcitriol": "calcitriol_LM",
    "calcipotriol": "calcipotriol_LM"
}


# ============================================================
# 3. VALIDATE AND BUILD HEATMAP MATRICES
# ============================================================

def build_matrix(data, column_map):

    required = [
        "cell_id",
        "hallmark",
        *column_map.values()
    ]

    missing = [
        col for col in required
        if col not in data.columns
    ]

    if missing:
        raise ValueError(
            f"Missing columns: {missing}"
        )

    matrix = np.empty((12, 3))
    records = []

    row = 0

    for cell in contexts:

        for hallmark in hallmarks:

            subset = data.loc[
                (data["cell_id"] == cell)
                &
                (data["hallmark"] == hallmark)
            ]

            if len(subset) != 1:

                raise ValueError(
                    f"Expected one row for "
                    f"{cell} / {hallmark}"
                )

            for j, treatment in enumerate(treatments):

                value = float(
                    subset[column_map[treatment]].iloc[0]
                )

                if not np.isfinite(value):
                    raise ValueError(
                        f"Invalid value for "
                        f"{cell}, {hallmark}, {treatment}"
                    )

                matrix[row, j] = value

                records.append({
                    "cell_id": cell,
                    "hallmark": hallmark,
                    "tratamiento": treatment,
                    "score": value
                })

            row += 1

    return matrix, pd.DataFrame(records)


matrix_full, full_long = build_matrix(
    full,
    {t: t for t in treatments}
)

matrix_landmark, landmark_long = build_matrix(
    landmark,
    landmark_columns
)


# ============================================================
# 4. VALIDATE ORIGINAL SUMMARY DATA
# ============================================================

fraction_col = "fraccion_hallmarks_negativos_3de3"

def summary_values(data, name):

    required = ["cell_id", fraction_col]

    missing = [
        c for c in required
        if c not in data.columns
    ]

    if missing:
        raise ValueError(
            f"Missing columns in {name}: {missing}"
        )

    if data["cell_id"].duplicated().any():
        raise ValueError(
            f"Duplicated cell IDs in {name}"
        )

    indexed = data.set_index("cell_id")

    if not set(contexts).issubset(indexed.index):
        raise ValueError(
            f"Missing contexts in {name}"
        )

    values = pd.to_numeric(
        indexed.loc[contexts, fraction_col],
        errors="raise"
    ).to_numpy(dtype=float)

    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or (values > 1).any()
    ):
        raise ValueError(
            f"Invalid fractions in {name}"
        )

    return values


values_full = summary_values(
    summary_full, "Full"
)

values_landmark = summary_values(
    summary_landmark, "Landmark"
)


# ============================================================
# 5. COMMON COLOUR SCALE
# ============================================================

max_abs = max(
    np.abs(matrix_full).max(),
    np.abs(matrix_landmark).max()
)

limit = np.ceil(max_abs * 10) / 10

norm = TwoSlopeNorm(
    vmin=-limit,
    vcenter=0,
    vmax=limit
)

cmap = plt.get_cmap("RdBu_r")


# ============================================================
# 6. PUBLICATION STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white"
})


# ============================================================
# 7. CREATE A4 LANDSCAPE FIGURE
# ============================================================

fig = plt.figure(
    figsize=(29.7 / 2.54, 21 / 2.54),
    facecolor="white"
)

# Upper heatmaps: leave room for long row labels
ax1 = fig.add_axes([
    0.155, 0.445, 0.285, 0.425
])

ax2 = fig.add_axes([
    0.625, 0.445, 0.285, 0.425
])

# Shared colour bar
cax = fig.add_axes([
    0.935, 0.445, 0.014, 0.425
])

# Lower summary panel
ax3 = fig.add_axes([
    0.135, 0.125, 0.77, 0.215
])


# ============================================================
# 8. PANEL LETTERS
# ============================================================

fig.text(
    0.025, 0.975,
    "A",
    fontsize=18,
    fontweight="bold",
    ha="left",
    va="top"
)

fig.text(
    0.495, 0.975,
    "B",
    fontsize=18,
    fontweight="bold",
    ha="left",
    va="top"
)

fig.text(
    0.025, 0.385,
    "C",
    fontsize=18,
    fontweight="bold",
    ha="left",
    va="top"
)


# ============================================================
# 9. TITLES
# ============================================================

fig.text(
    0.2975, 0.925,
    "Full L1000 gene space",
    ha="center",
    va="center",
    fontsize=14,
    fontweight="bold"
)

fig.text(
    0.7675, 0.925,
    "Landmark genes",
    ha="center",
    va="center",
    fontsize=14,
    fontweight="bold"
)

fig.text(
    0.52, 0.375,
    "Context-level convergence of the core programmes",
    ha="center",
    va="center",
    fontsize=13,
    fontweight="bold"
)


# ============================================================
# 10. HEATMAPS
# ============================================================

row_labels = [
    f"{cell} | {hallmark_labels[hallmark]}"
    for cell in contexts
    for hallmark in hallmarks
]


def draw_heatmap(ax, matrix):

    image = ax.imshow(
        matrix,
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest"
    )

    ax.set_xticks(range(3))

    ax.set_xticklabels(
        [treatment_labels[t] for t in treatments],
        fontsize=10.5
    )

    ax.set_yticks(range(12))

    ax.set_yticklabels(
        row_labels,
        fontsize=9.5
    )

    ax.tick_params(
        axis="both",
        length=0,
        pad=5
    )

    # White cell boundaries
    ax.set_xticks(
        np.arange(-0.5, 3, 1),
        minor=True
    )

    ax.set_yticks(
        np.arange(-0.5, 12, 1),
        minor=True
    )

    ax.grid(
        which="minor",
        color="white",
        linewidth=0.9
    )

    ax.tick_params(
        which="minor",
        bottom=False,
        left=False
    )

    # Stronger separation between cellular contexts
    for boundary in [2.5, 5.5, 8.5]:

        ax.axhline(
            boundary,
            color="#444444",
            linewidth=1.2
        )

    for i in range(12):

        for j in range(3):

            value = matrix[i, j]

            rgba = cmap(norm(value))

            luminance = (
                0.2126 * rgba[0]
                + 0.7152 * rgba[1]
                + 0.0722 * rgba[2]
            )

            text_colour = (
                "#111111"
                if luminance > 0.52
                else "white"
            )

            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=9.3,
                fontweight="medium",
                color=text_colour
            )

    return image


im1 = draw_heatmap(
    ax1,
    matrix_full
)

im2 = draw_heatmap(
    ax2,
    matrix_landmark
)


# ============================================================
# 11. COMMON COLOUR BAR
# ============================================================

cbar = fig.colorbar(
    im2,
    cax=cax
)

cbar.set_label(
    "Mean Hallmark programme score",
    fontsize=10.5,
    labelpad=12
)

cbar.ax.tick_params(
    labelsize=9.5
)


# ============================================================
# 12. PANEL C — CONTEXT CONVERGENCE
# ============================================================

x = np.arange(len(contexts))
width = 0.34

bars_full = ax3.bar(
    x - width / 2,
    values_full,
    width,
    color="#0072B2",
    edgecolor="black",
    linewidth=0.5,
    label="Full",
    zorder=3
)

bars_landmark = ax3.bar(
    x + width / 2,
    values_landmark,
    width,
    color="#E69F00",
    edgecolor="black",
    linewidth=0.5,
    label="Landmark",
    zorder=3
)

ax3.set_xticks(x)

ax3.set_xticklabels(
    contexts,
    fontsize=11
)

ax3.set_ylim(0, 1.16)

ax3.set_ylabel(
    "Fraction of core Hallmarks\n"
    "negative in all 3 treatments",
    fontsize=10.5
)

ax3.legend(
    loc="upper left",
    frameon=False,
    fontsize=10,
    ncol=2
)

ax3.grid(
    axis="y",
    color="#E6E6E6",
    linewidth=0.6,
    zorder=0
)

ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)


# ============================================================
# 13. FRACTION LABELS
# ============================================================

for bars, values in [
    (bars_full, values_full),
    (bars_landmark, values_landmark)
]:

    for bar, value in zip(bars, values):

        count = round(value * 3)

        if not np.isclose(
            value, count / 3,
            atol=1e-6
        ):
            raise ValueError(
                "Fraction is not a multiple of 1/3."
            )

        ax3.annotate(
            f"{count}/3",
            (
                bar.get_x()
                + bar.get_width() / 2,
                bar.get_height()
            ),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold"
        )


# ============================================================
# 14. METHODOLOGICAL NOTE
# ============================================================

fig.text(
    0.52, 0.045,
    "Context-level values are descriptive because cell line "
    "and experimental block are confounded in the matched design.",
    ha="center",
    va="center",
    fontsize=10,
    color="#333333"
)


# ============================================================
# 15. EXPORT — OVERWRITE PREVIOUS FILES
# ============================================================

OUTPUT.mkdir(
    parents=True,
    exist_ok=True
)

pdf_path = (
    OUTPUT / "Figure_5_context_dependence.pdf"
)

png_path = (
    OUTPUT / "Figure_5_context_dependence.png"
)

csv_path = (
    OUTPUT / "Figure_5_context_values.csv"
)

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

audit_table = pd.concat(
    [
        full_long.assign(
            analysis_space="Full"
        ),
        landmark_long.assign(
            analysis_space="Landmark"
        )
    ],
    ignore_index=True
)

audit_table.to_csv(
    csv_path,
    index=False
)

plt.close(fig)

print("\n" + "=" * 75)
print("FIGURE 5 GENERATED SUCCESSFULLY")
print("=" * 75)
print("Size: A4 landscape — 29.7 x 21 cm")
print("PDF:", pdf_path)
print("PNG:", png_path)
print("Audit CSV:", csv_path)
