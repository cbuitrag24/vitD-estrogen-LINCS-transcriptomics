
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.colors import TwoSlopeNorm

# ============================================================
# FIGURE 4 — HALLMARK CORE CONVERGENCE
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

FILENAME = "GSEA_FULL_vs_LANDMARK_comparacion.csv"

print("=" * 75)
print("FIGURE 4 — HALLMARK CORE CONVERGENCE")
print("=" * 75)


# ============================================================
# 1. LOCATE ORIGINAL RESULTS
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


data_path = locate_file(FILENAME)
print("Source:", data_path)

df = pd.read_csv(data_path)


# ============================================================
# 2. SELECT ORIGINAL HALLMARK RESULTS
# ============================================================

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

required = [
    "tratamiento",
    "hallmark",
    "NES_full",
    "q_full",
    "NES_landmark",
    "q_landmark",
    "robustez"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(f"Missing columns: {missing}")

core = df.loc[
    df["hallmark"].isin(hallmarks)
    &
    df["tratamiento"].isin(treatments)
].copy()

expected = {
    (hallmark, treatment)
    for hallmark in hallmarks
    for treatment in treatments
}

actual = list(zip(
    core["hallmark"],
    core["tratamiento"]
))

if (
    len(actual) != 9
    or len(set(actual)) != 9
    or set(actual) != expected
):
    raise ValueError(
        "Expected exactly nine unique Hallmark-treatment combinations."
    )

numeric = [
    "NES_full",
    "q_full",
    "NES_landmark",
    "q_landmark"
]

for col in numeric:
    core[col] = pd.to_numeric(
        core[col], errors="raise"
    )

if not np.isfinite(core[numeric].to_numpy()).all():
    raise ValueError("Non-finite NES or q values.")

if (
    (core[["q_full", "q_landmark"]] < 0).any().any()
    or
    (core[["q_full", "q_landmark"]] > 1).any().any()
):
    raise ValueError("Adjusted P values outside [0, 1].")


# ============================================================
# 3. CONSTRUCT MATRICES
# ============================================================

nes_full = np.empty((3, 3))
nes_landmark = np.empty((3, 3))

q_full = np.empty((3, 3))
q_landmark = np.empty((3, 3))

for i, hallmark in enumerate(hallmarks):
    for j, treatment in enumerate(treatments):

        row = core.loc[
            (core["hallmark"] == hallmark)
            &
            (core["tratamiento"] == treatment)
        ].iloc[0]

        nes_full[i, j] = row["NES_full"]
        nes_landmark[i, j] = row["NES_landmark"]

        q_full[i, j] = row["q_full"]
        q_landmark[i, j] = row["q_landmark"]


# ============================================================
# 4. COMMON COLOUR SCALE
# ============================================================

max_abs = max(
    np.abs(nes_full).max(),
    np.abs(nes_landmark).max()
)

limit = np.ceil(max_abs * 10) / 10

norm = TwoSlopeNorm(
    vmin=-limit,
    vcenter=0,
    vmax=limit
)

cmap = plt.get_cmap("RdBu_r")


# ============================================================
# 5. PUBLICATION STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 15,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


# ============================================================
# 6. A4 LANDSCAPE FIGURE
# ============================================================

fig = plt.figure(
    figsize=(29.7 / 2.54, 21 / 2.54),
    facecolor="white"
)

# Large panels with room for labels and colour bar
ax1 = fig.add_axes([0.105, 0.23, 0.34, 0.58])
ax2 = fig.add_axes([0.565, 0.23, 0.34, 0.58])

cax = fig.add_axes([0.93, 0.23, 0.018, 0.58])

# Panel letters, outside heatmaps
fig.text(
    0.025, 0.96, "A",
    fontsize=20,
    fontweight="bold",
    ha="left",
    va="top"
)

fig.text(
    0.49, 0.96, "B",
    fontsize=20,
    fontweight="bold",
    ha="left",
    va="top"
)

# Titles in separate row
fig.text(
    0.275, 0.895,
    "Full L1000 gene space",
    fontsize=15,
    fontweight="bold",
    ha="center",
    va="center"
)

fig.text(
    0.735, 0.895,
    "Landmark genes",
    fontsize=15,
    fontweight="bold",
    ha="center",
    va="center"
)


# ============================================================
# 7. SIGNIFICANCE LABELS
# ============================================================

def stars(q):
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


# ============================================================
# 8. HEATMAP FUNCTION
# ============================================================

def plot_heatmap(ax, nes, q_values):

    image = ax.imshow(
        nes,
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest"
    )

    ax.set_xticks(range(3))
    ax.set_xticklabels(
        [treatment_labels[t] for t in treatments],
        fontsize=12
    )

    ax.set_yticks(range(3))
    ax.set_yticklabels(
        [hallmark_labels[h] for h in hallmarks],
        fontsize=12
    )

    ax.tick_params(
        axis="both",
        length=0,
        pad=10
    )

    # Thin white boundaries between cells
    ax.set_xticks(
        np.arange(-0.5, 3, 1),
        minor=True
    )

    ax.set_yticks(
        np.arange(-0.5, 3, 1),
        minor=True
    )

    ax.grid(
        which="minor",
        color="white",
        linewidth=2
    )

    ax.tick_params(
        which="minor",
        bottom=False,
        left=False
    )

    for i in range(3):
        for j in range(3):

            value = nes[i, j]
            q = q_values[i, j]

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
                f"NES {value:.2f}{stars(q)}",
                ha="center",
                va="center",
                fontsize=12,
                fontweight="bold",
                color=text_colour
            )

    return image


# ============================================================
# 9. PLOT BOTH PANELS
# ============================================================

im1 = plot_heatmap(
    ax1,
    nes_full,
    q_full
)

im2 = plot_heatmap(
    ax2,
    nes_landmark,
    q_landmark
)


# ============================================================
# 10. COMMON COLOUR BAR
# ============================================================

cbar = fig.colorbar(
    im2,
    cax=cax
)

cbar.set_label(
    "Normalised enrichment score (NES)",
    fontsize=12,
    labelpad=13
)

cbar.ax.tick_params(
    labelsize=11
)


# ============================================================
# 11. STATISTICAL NOTE
# ============================================================

fig.text(
    0.50,
    0.095,
    "Negative NES indicates negative enrichment. "
    "* BH-adjusted empirical P < 0.05; "
    "** P < 0.01; "
    "*** P < 0.001.",
    ha="center",
    va="center",
    fontsize=11,
    color="#333333"
)


# ============================================================
# 12. EXPORT — OVERWRITE PREVIOUS FILES
# ============================================================

OUTPUT.mkdir(
    parents=True,
    exist_ok=True
)

pdf_path = (
    OUTPUT / "Figure_4_Hallmark_core_convergence.pdf"
)

png_path = (
    OUTPUT / "Figure_4_Hallmark_core_convergence.png"
)

csv_path = (
    OUTPUT / "Figure_4_Hallmark_core_values.csv"
)

# Preserve exact A4 size
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

audit_table = core[
    [
        "tratamiento",
        "hallmark",
        "NES_full",
        "q_full",
        "NES_landmark",
        "q_landmark",
        "robustez"
    ]
].copy()

audit_table.to_csv(
    csv_path,
    index=False
)

plt.close(fig)

print("\nFIGURE 4 GENERATED")
print("Size: A4 landscape — 29.7 x 21 cm")
print("PDF:", pdf_path)
print("PNG:", png_path)
print("Audit CSV:", csv_path)
