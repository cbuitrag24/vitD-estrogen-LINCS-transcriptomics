
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ============================================================
# FIGURE 2 — CORRECTED PUBLICATION LAYOUT
# Exact size: 14.5 x 9 cm
# PDF vector + PNG 600 dpi
# ============================================================

SOURCE = Path(
    r"C:\Users\petri\LINCS_Estradiol_VitaminaD"
)

REPO = Path(
    r"C:\Users\petri\vitD-estrogen-LINCS-transcriptomics"
)

OUTPUT = REPO / "results" / "figures"

# ------------------------------------------------------------
# 1. Original PERMANOVA results
# ------------------------------------------------------------

permanova = pd.DataFrame({
    "Gene space": [
        "Full", "Full", "Landmark", "Landmark"
    ],
    "Source": [
        "Treatment", "Context / block",
        "Treatment", "Context / block"
    ],
    "R2": [
        0.1380151871065361,
        0.4252446816912371,
        0.14877538546041968,
        0.42440058382449297
    ]
})

p_values = [0.5740740740740741, 0.2916666666666667]

# ------------------------------------------------------------
# 2. Locate original Pearson similarity data
# ------------------------------------------------------------

filename = "ESTADISTICA_similitud_por_contexto.csv"

matches = []
for folder in (SOURCE, REPO):
    if folder.is_dir():
        matches.extend(folder.rglob(filename))

matches = list({
    str(p.resolve()): p
    for p in matches if p.is_file()
}.values())

if len(matches) != 1:
    raise FileNotFoundError(
        f"Expected one original similarity CSV; "
        f"found {len(matches)}:\n" +
        "\n".join(str(p) for p in matches)
    )

sim = pd.read_csv(matches[0])

required = [
    "analysis", "metric", "cell_id",
    "vitd_vitd", "e2_calcitriol", "e2_calcipotriol"
]

missing = [c for c in required if c not in sim.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

sim = sim.loc[
    sim["analysis"].astype(str).str.lower()
    .str.contains("full", na=False)
    &
    sim["metric"].astype(str).str.lower()
    .str.contains("pearson", na=False)
].copy()

comparison_codes = [
    "vitd_vitd",
    "e2_calcitriol",
    "e2_calcipotriol"
]

sim_long = sim.melt(
    id_vars=["cell_id"],
    value_vars=comparison_codes,
    var_name="comparison",
    value_name="pearson"
)

sim_long["cell_id"] = (
    sim_long["cell_id"].astype(str).str.strip().str.upper()
)

sim_long["pearson"] = pd.to_numeric(
    sim_long["pearson"], errors="raise"
)

cells = ["A549", "HA1E", "PC3", "VCAP"]

expected = {
    (cell, comparison)
    for cell in cells
    for comparison in comparison_codes
}

actual = list(zip(
    sim_long["cell_id"], sim_long["comparison"]
))

if (
    len(actual) != 12
    or set(actual) != expected
    or len(set(actual)) != 12
):
    raise ValueError(
        "Expected exactly one observation for each "
        "cell line and comparison."
    )

if (
    not np.isfinite(sim_long["pearson"]).all()
    or (sim_long["pearson"].abs() > 1).any()
):
    raise ValueError("Invalid Pearson correlations")

means = [
    sim_long.loc[
        sim_long["comparison"] == code, "pearson"
    ].mean()
    for code in comparison_codes
]

print("Pearson means:", [round(x, 4) for x in means])

# ------------------------------------------------------------
# 3. Figure style and fixed layout
# ------------------------------------------------------------

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

fig = plt.figure(
    figsize=(14.5 / 2.54, 9 / 2.54),
    facecolor="white"
)

# Axes are positioned explicitly.
# No constrained_layout or tight_layout.
ax1 = fig.add_axes([0.085, 0.28, 0.34, 0.46])
ax2 = fig.add_axes([0.545, 0.28, 0.415, 0.46])

# Panel labels and titles occupy separate rows.
fig.text(
    0.035, 0.955, "A",
    fontsize=12, fontweight="bold",
    ha="left", va="top"
)

fig.text(
    0.50, 0.955, "B",
    fontsize=12, fontweight="bold",
    ha="left", va="top"
)

fig.text(
    0.255, 0.92,
    "Global variance partitioning",
    fontsize=9.5, fontweight="bold",
    ha="center", va="center"
)

fig.text(
    0.752, 0.92,
    "Pairwise transcriptional similarity\nacross contexts",
    fontsize=9.5, fontweight="bold",
    ha="center", va="top"
)

# ============================================================
# PANEL A — Variance partitioning
# ============================================================

x = np.arange(2)
width = 0.30

r2_t = [0.1380151871065361, 0.14877538546041968]
r2_c = [0.4252446816912371, 0.42440058382449297]

blue = "#0072B2"
orange = "#E69F00"

bars_t = ax1.bar(
    x - width / 2, r2_t, width,
    color=blue, edgecolor="black", linewidth=0.35
)

bars_c = ax1.bar(
    x + width / 2, r2_c, width,
    color=orange, edgecolor="black", linewidth=0.35
)

ax1.set_xticks(x, ["Full", "Landmark"])
ax1.set_ylabel("Variance explained (R²)")
ax1.set_ylim(0, 0.55)

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

for bars in (bars_t, bars_c):
    for bar in bars:
        ax1.annotate(
            f"{bar.get_height():.2f}",
            (
                bar.get_x() + bar.get_width() / 2,
                bar.get_height()
            ),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center", va="bottom",
            fontsize=8
        )

# p-values placed below the treatment bars,
# outside the plot, in a separate line.
for i, (bar, p) in enumerate(zip(bars_t, p_values)):
    ax1.annotate(
        f"p = {p:.3f}",
        (
            bar.get_x() + bar.get_width() / 2,
            0
        ),
        xytext=(0, -27),
        textcoords="offset points",
        ha="center", va="top",
        fontsize=7.4,
        annotation_clip=False
    )

# Treatment legend, clearly separated from title.
fig.legend(
    handles=[
        Patch(facecolor=blue, label="Treatment"),
        Patch(facecolor=orange, label="Context / block")
    ],
    loc="center",
    bbox_to_anchor=(0.255, 0.82),
    ncol=2,
    frameon=False,
    fontsize=7.7,
    handlelength=1.0,
    columnspacing=1.0
)

# ============================================================
# PANEL B — Pearson correlations
# ============================================================

markers = {
    "A549": "o",
    "HA1E": "s",
    "PC3": "^",
    "VCAP": "D"
}

colors = {
    "A549": "#0072B2",
    "HA1E": "#E69F00",
    "PC3": "#009E73",
    "VCAP": "#CC79A7"
}

offsets = {
    "A549": -0.18,
    "HA1E": -0.06,
    "PC3": 0.06,
    "VCAP": 0.18
}

for i, code in enumerate(comparison_codes):
    subset = sim_long.loc[
        sim_long["comparison"] == code
    ]

    for cell in cells:
        row = subset.loc[
            subset["cell_id"] == cell
        ].iloc[0]

        ax2.scatter(
            i + offsets[cell],
            row["pearson"],
            s=55,
            marker=markers[cell],
            color=colors[cell],
            edgecolors="black",
            linewidths=0.5,
            zorder=3
        )

    ax2.plot(
        [i - 0.25, i + 0.25],
        [means[i], means[i]],
        color="black",
        linewidth=1.8,
        zorder=2
    )

    ax2.annotate(
        f"Mean = {means[i]:.2f}",
        (i, means[i]),
        xytext=(0, -18) if i == 1 else (0, 7),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=7.3,
        fontweight="bold",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.85,
            pad=1
        )
    )

ax2.axhline(
    0, color="#999999",
    linewidth=0.7, zorder=1
)

ax2.set_xlim(-0.42, 2.42)
ax2.set_ylim(-0.16, 0.63)

ax2.set_xticks([0, 1, 2])
ax2.set_xticklabels([
    "Calcitriol vs\nCalcipotriol",
    "Estradiol vs\nCalcitriol",
    "Estradiol vs\nCalcipotriol"
], fontsize=7.3)

ax2.set_ylabel("Pearson correlation")

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# ------------------------------------------------------------
# 4. Separate notes and legend
# ------------------------------------------------------------

fig.text(
    0.255, 0.105,
    "Treatment: exact within-block permutations\n"
    "Context/block R²: descriptive",
    ha="center",
    va="center",
    fontsize=7.2,
    color="#444444"
)

handles = [
    Line2D(
        [0], [0],
        marker=markers[cell],
        linestyle="None",
        markerfacecolor=colors[cell],
        markeredgecolor="black",
        markeredgewidth=0.5,
        markersize=6,
        label=cell
    )
    for cell in cells
]

fig.text(
    0.752, 0.165,
    "Cell line",
    ha="center",
    va="center",
    fontsize=8,
    fontweight="bold"
)

fig.legend(
    handles=handles,
    loc="center",
    bbox_to_anchor=(0.752, 0.105),
    ncol=4,
    frameon=False,
    fontsize=7.5,
    columnspacing=0.85,
    handletextpad=0.3
)

# ------------------------------------------------------------
# 5. Save and overwrite
# ------------------------------------------------------------

OUTPUT.mkdir(parents=True, exist_ok=True)

pdf_path = OUTPUT / "Figure_2_global_variance_similarity.pdf"
png_path = OUTPUT / "Figure_2_global_variance_similarity.png"
csv_path = OUTPUT / "Figure_2_PERMANOVA_values.csv"

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

permanova.to_csv(csv_path, index=False)

plt.close(fig)

print("\nFIGURE 2 GENERATED")
print("Size: 14.5 x 9 cm")
print("PDF:", pdf_path)
print("PNG:", png_path)
print("CSV:", csv_path)
