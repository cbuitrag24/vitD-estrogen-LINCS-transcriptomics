
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from sklearn.decomposition import PCA
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

# ============================================================
# FIGURE 1 — FINAL LAYOUT REVISION
# Size: 14.5 x 9 cm
# PDF vector + PNG 600 dpi
# ============================================================

SOURCE = Path(
    r"C:\Users\petri\LINCS_Estradiol_VitaminaD"
)

REPO = Path(
    r"C:\Users\petri\vitD-estrogen-LINCS-transcriptomics"
)

OUTPUT = REPO / "results" / "figures"

NPZ = SOURCE / "dataset_primario_expresion_12x12328.npz"

META_NAMES = [
    "dataset_primario_metadata.csv",
    "metadata_primario.csv",
    "diseno_primario.csv",
    "firmas_primarias.csv",
]

# ------------------------------------------------------------
# 1. Load original data
# ------------------------------------------------------------

if not NPZ.is_file():
    raise FileNotFoundError(f"Missing file: {NPZ}")

meta_path = next(
    (
        SOURCE / name
        for name in META_NAMES
        if (SOURCE / name).is_file()
    ),
    None
)

if meta_path is None:
    raise FileNotFoundError(
        "Primary metadata CSV not found in source folder."
    )

with np.load(NPZ) as archive:
    X = archive["expresion"]

meta = pd.read_csv(meta_path)

if X.shape != (12, 12328):
    raise ValueError(f"Unexpected matrix shape: {X.shape}")

if len(meta) != 12:
    raise ValueError("Metadata must contain 12 rows")

if not np.isfinite(X).all():
    raise ValueError("Expression matrix has invalid values")

def find_column(names):
    for name in names:
        for col in meta.columns:
            if col.lower() == name.lower():
                return col
    return None

cell_col = find_column(
    ["cell_id", "cell", "cell_line", "celula"]
)

treat_col = find_column(
    ["pert_iname", "treatment", "tratamiento", "compound"]
)

if cell_col is None or treat_col is None:
    raise ValueError("Missing cell or treatment column")

# ------------------------------------------------------------
# 2. PCA — original analytical method
# ------------------------------------------------------------

pca = PCA(n_components=2)
coords = pca.fit_transform(X)
variance = pca.explained_variance_ratio_ * 100

print(f"PC1 = {variance[0]:.2f}%")
print(f"PC2 = {variance[1]:.2f}%")

if (
    abs(variance[0] - 33.94) >= 0.15
    or abs(variance[1] - 15.23) >= 0.15
):
    raise ValueError(
        "PCA variance differs from manuscript values. "
        "No figure will be overwritten."
    )

def treatment_name(value):
    s = str(value).lower()
    if "estradiol" in s:
        return "Estradiol"
    if "calcipotriol" in s:
        return "Calcipotriol"
    if "calcitriol" in s:
        return "Calcitriol"
    return str(value)

df = pd.DataFrame({
    "PC1": coords[:, 0],
    "PC2": coords[:, 1],
    "cell_id": meta[cell_col].astype(str).str.strip(),
    "treatment": meta[treat_col].astype(str),
})

df["treatment_label"] = df["treatment"].map(
    treatment_name
)

cell_numbers = {
    "A549": 1,
    "HA1E": 2,
    "PC3": 3,
    "VCAP": 4,
}

df["cell_number"] = (
    df["cell_id"].str.upper().map(cell_numbers)
)

if df["cell_number"].isna().any():
    raise ValueError(
        "Unexpected cell-line name in metadata"
    )

if set(df["treatment_label"]) != {
    "Estradiol", "Calcitriol", "Calcipotriol"
}:
    raise ValueError("Unexpected treatments")

if not (
    df["treatment_label"].value_counts() == 4
).all():
    raise ValueError(
        "Expected four samples per treatment"
    )

if not (
    df.groupby(["treatment_label", "cell_number"])
    .size()
    .eq(1)
    .all()
    and len(df.groupby(
        ["treatment_label", "cell_number"]
    )) == 12
):
    raise ValueError(
        "Expected one sample per cell/treatment combination"
    )

# ------------------------------------------------------------
# 3. Figure setup
# ------------------------------------------------------------

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.7,
})

fig = plt.figure(
    figsize=(14.5 / 2.54, 9 / 2.54),
    facecolor="white"
)

gs = fig.add_gridspec(
    1, 2,
    width_ratios=[1.05, 1.45],
    left=0.055,
    right=0.96,
    top=0.78,
    bottom=0.30,
    wspace=0.38
)

# ============================================================
# PANEL A
# ============================================================

ax1 = fig.add_subplot(gs[0, 0])
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax1.axis("off")

ax1.text(
    0.00, 1.37, "A",
    transform=ax1.transAxes,
    fontsize=12,
    fontweight="bold",
    va="top"
)

ax1.text(
    0.50, 1.28,
    "Matched LINCS L1000\nstudy design",
    transform=ax1.transAxes,
    ha="center",
    va="top",
    fontsize=9.4,
    fontweight="bold"
)

def box(y, text):
    width = 0.98
    height = 0.16

    patch = FancyBboxPatch(
        (0.01, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        transform=ax1.transAxes,
        facecolor="#F7F8FA",
        edgecolor="#555555",
        linewidth=0.75,
        clip_on=False
    )
    ax1.add_patch(patch)

    ax1.text(
        0.5, y, text,
        transform=ax1.transAxes,
        ha="center",
        va="center",
        fontsize=7.3,
        linespacing=1.18
    )

def arrow(y1, y2):
    ax1.annotate(
        "",
        xy=(0.5, y2),
        xytext=(0.5, y1),
        xycoords="axes fraction",
        textcoords="axes fraction",
        arrowprops=dict(
            arrowstyle="-|>",
            color="#555555",
            lw=0.85,
            mutation_scale=9
        )
    )

box(
    0.87,
    "LINCS L1000 signatures\n10 µM · 24 h"
)
arrow(0.78, 0.68)

box(
    0.59,
    "4 matched experimental blocks\nA549 · HA1E · PC3 · VCAP"
)
arrow(0.50, 0.40)

box(
    0.31,
    "Estradiol · Calcitriol\nCalcipotriol"
)
arrow(0.22, 0.12)

box(
    0.03,
    "12 Level 5 MODZ signatures\n4 blocks × 3 treatments"
)

# ============================================================
# PANEL B — PCA
# ============================================================

ax2 = fig.add_subplot(gs[0, 1])

ax2.text(
    -0.12, 1.37, "B",
    transform=ax2.transAxes,
    fontsize=12,
    fontweight="bold",
    va="top"
)

ax2.text(
    0.5, 1.28,
    "Global transcriptional\nstructure",
    transform=ax2.transAxes,
    ha="center",
    va="top",
    fontsize=9.4,
    fontweight="bold"
)

styles = {
    "Estradiol": ("o", "#0072B2"),
    "Calcitriol": ("s", "#E69F00"),
    "Calcipotriol": ("^", "#009E73"),
}

# Use treatment-specific marker and colour.
# Numbers 1–4 identify cell lines.

for treatment, (marker, color) in styles.items():

    subset = df[
        df["treatment_label"] == treatment
    ]

    ax2.scatter(
        subset["PC1"],
        subset["PC2"],
        s=110,
        marker=marker,
        c=color,
        edgecolors="black",
        linewidths=0.6,
        zorder=3
    )

    # Each number is drawn inside its corresponding symbol.
    for _, row in subset.iterrows():
        ax2.text(
            row["PC1"],
            row["PC2"],
            str(int(row["cell_number"])),
            fontsize=6.2,
            fontweight="bold",
            color="black",
            ha="center",
            va="center",
            zorder=5
        )

ax2.axhline(
    0, color="#C8C8C8", lw=0.7, zorder=1
)

ax2.axvline(
    0, color="#C8C8C8", lw=0.7, zorder=1
)

ax2.set_xlabel(
    f"PC1 ({variance[0]:.1f}% variance)",
    labelpad=4
)

ax2.set_ylabel(
    f"PC2 ({variance[1]:.1f}% variance)",
    labelpad=4
)

ax2.tick_params(
    direction="out",
    length=3,
    width=0.7
)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Ensure all markers remain inside axes.
for key in ("PC1", "PC2"):
    low = df[key].min()
    high = df[key].max()
    padding = max((high - low) * 0.13, 1.0)

    if key == "PC1":
        ax2.set_xlim(low - padding, high + padding)
    else:
        ax2.set_ylim(low - padding, high + padding)

# ------------------------------------------------------------
# 4. Legends within the figure canvas
# ------------------------------------------------------------

handles = [
    Line2D(
        [0], [0],
        marker=styles[t][0],
        linestyle="None",
        markerfacecolor=styles[t][1],
        markeredgecolor="black",
        markeredgewidth=0.5,
        markersize=6.5,
        label=t
    )
    for t in styles
]

# Treatment legend: below the PCA axis.
legend = ax2.legend(
    handles=handles,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.32),
    ncol=3,
    frameon=False,
    fontsize=7.1,
    handletextpad=0.35,
    columnspacing=0.85,
    borderaxespad=0
)

# Cell-line key: second row below treatment legend.
ax2.text(
    0.5, -0.56,
    "Cell lines:  1 A549   ·   2 HA1E   ·   "
    "3 PC3   ·   4 VCAP",
    transform=ax2.transAxes,
    ha="center",
    va="top",
    fontsize=7.2,
    clip_on=False
)

# ------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------

OUTPUT.mkdir(parents=True, exist_ok=True)

pdf_path = OUTPUT / "Figure_1_study_design_PCA.pdf"
png_path = OUTPUT / "Figure_1_study_design_PCA.png"
csv_path = OUTPUT / "Figure_1_PCA_coordinates.csv"

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

df.to_csv(csv_path, index=False)

plt.close(fig)

print("=" * 65)
print("FIGURE 1 GENERATED")
print("Dimensions: 14.5 x 9 cm")
print("PDF:", pdf_path)
print("PNG:", png_path)
print("CSV:", csv_path)
print("=" * 65)
