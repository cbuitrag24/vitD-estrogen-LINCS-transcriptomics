
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.lines import Line2D

# ============================================================
# FIGURE 6 — TEMPORAL SENSITIVITY
# A4 LANDSCAPE — PUBLICATION VERSION
# ============================================================

# Resolve all paths relative to this script, independent of the user's machine.
REPO = Path(__file__).resolve().parent.parent
SOURCE = REPO  # Scripts 50 and 51 write these intermediate CSVs to the repository root.

OUTPUT = REPO / "results" / "figures"
OUTPUT.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("FIGURE 6 — IMPROVED PUBLICATION LAYOUT")
print("=" * 70)


# ============================================================
# 1. LOAD DATA
# ============================================================

def load_csv(name):
    path = SOURCE / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


sim_full = load_csv("E2_temporal_similitud_global.csv")
sim_lm = load_csv("E2_temporal_similitud_landmark.csv")
hall_full = load_csv("E2_temporal_5Hallmarks_por_bloque.csv")
hall_lm = load_csv(
    "E2_temporal_5Hallmarks_landmark_por_bloque.csv"
)

keys = ["cell_id", "bloque"]

hallmarks = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1"
]

labels = [
    "E2F TARGETS",
    "G2M CHECKPOINT",
    "MYC TARGETS V1"
]


# ============================================================
# 2. MATCH BLOCKS
# ============================================================

sim = sim_full[
    keys + ["pearson_6h_24h"]
].merge(
    sim_lm[keys + ["pearson"]],
    on=keys,
    how="outer",
    validate="one_to_one",
    indicator=True
)

if not (sim["_merge"] == "both").all():
    raise ValueError("Similarity blocks do not match.")

sim = sim.drop(columns="_merge").rename(columns={
    "pearson_6h_24h": "Full",
    "pearson": "Landmark"
})

sim = sim.sort_values(keys).reset_index(drop=True)

if len(sim) != 7:
    raise ValueError("Expected seven matched blocks.")

if not np.isfinite(
    sim[["Full", "Landmark"]].to_numpy(dtype=float)
).all():
    raise ValueError("Non-finite similarity values.")


# ============================================================
# 3. PREPARE ORIGINAL HALLMARK VALUES
# ============================================================

def prepare(data, space):

    subset = data.loc[
        data["hallmark"].isin(hallmarks)
    ].copy()

    required = [
        "hallmark",
        "cell_id",
        "bloque",
        "score_6h",
        "score_24h",
        "delta_24h_minus_6h"
    ]

    if not set(required).issubset(subset.columns):
        raise ValueError(f"Missing columns: {space}")

    if subset.duplicated(
        ["hallmark", "cell_id", "bloque"]
    ).any():
        raise ValueError(f"Duplicated blocks: {space}")

    expected = {
        (h, c, b)
        for h in hallmarks
        for c, b in sim[keys].itertuples(
            index=False, name=None
        )
    }

    observed = set(
        subset[
            ["hallmark", "cell_id", "bloque"]
        ].itertuples(index=False, name=None)
    )

    if expected != observed:
        raise ValueError(
            f"Hallmark/block coverage mismatch: {space}"
        )

    if not np.allclose(
        subset["score_24h"] - subset["score_6h"],
        subset["delta_24h_minus_6h"],
        atol=1e-5,
        rtol=0
    ):
        raise ValueError(
            f"Inconsistent temporal differences: {space}"
        )

    subset["analysis_space"] = space
    return subset


hf = prepare(hall_full, "Full")
hl = prepare(hall_lm, "Landmark")

combined = pd.concat([hf, hl], ignore_index=True)


# ============================================================
# 4. TYPOGRAPHY
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 15,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "axes.linewidth": 0.9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

BLUE = "#0072B2"
ORANGE = "#E69F00"

fig = plt.figure(
    figsize=(29.7 / 2.54, 21 / 2.54),
    facecolor="white"
)

# Two spacious panels with room for labels
ax_a = fig.add_axes([0.12, 0.61, 0.82, 0.28])
ax_b = fig.add_axes([0.12, 0.18, 0.82, 0.28])


# ============================================================
# 5. PANEL A — GLOBAL SIMILARITY
# ============================================================

x = np.arange(len(sim))
width = 0.35

ax_a.bar(
    x - width / 2,
    sim["Full"],
    width,
    color=BLUE,
    label="Full",
    zorder=3
)

ax_a.bar(
    x + width / 2,
    sim["Landmark"],
    width,
    color=ORANGE,
    label="Landmark",
    zorder=3
)

ax_a.axhline(0, color="black", linewidth=0.9)

# Two-line labels: easier to read than rotated labels
block_labels = [
    f"{cell}\n{block}"
    for cell, block in sim[keys].itertuples(
        index=False, name=None
    )
]

ax_a.set_xticks(x)
ax_a.set_xticklabels(
    block_labels,
    fontsize=11,
    rotation=0,
    ha="center"
)

ax_a.tick_params(axis="x", pad=8)

ax_a.set_ylabel(
    "Pearson correlation\n(6 h vs 24 h)",
    fontsize=12
)

ax_a.set_title(
    "Global transcriptional similarity",
    fontsize=15,
    fontweight="bold",
    pad=12
)

ax_a.legend(
    frameon=False,
    ncol=2,
    loc="upper left"
)

ax_a.grid(
    axis="y",
    color="#E5E5E5",
    linewidth=0.7,
    zorder=0
)

ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)

fig.text(
    0.035, 0.965,
    "A",
    fontsize=20,
    fontweight="bold",
    va="top"
)


# ============================================================
# 6. PANEL B — PROGRAMME-LEVEL TEMPORAL CHANGES
# ============================================================

positions = np.arange(3)

offsets = {
    "Full": -0.15,
    "Landmark": 0.15
}

colours = {
    "Full": BLUE,
    "Landmark": ORANGE
}

# Fixed offsets: no random scientific values or jitter
point_offsets = np.linspace(-0.07, 0.07, 7)

for i, hallmark in enumerate(hallmarks):

    for space in ["Full", "Landmark"]:

        subset = combined.loc[
            (combined["hallmark"] == hallmark)
            &
            (combined["analysis_space"] == space)
        ]

        subset = sim[keys].merge(
            subset,
            on=keys,
            how="left",
            validate="one_to_one"
        )

        values = subset[
            "delta_24h_minus_6h"
        ].to_numpy(dtype=float)

        centre = positions[i] + offsets[space]

        ax_b.scatter(
            centre + point_offsets,
            values,
            s=55,
            color=colours[space],
            alpha=0.80,
            edgecolor="white",
            linewidth=0.6,
            zorder=4
        )

        mean_value = np.mean(values)

        ax_b.plot(
            [centre - 0.11, centre + 0.11],
            [mean_value, mean_value],
            color=colours[space],
            linewidth=4,
            zorder=5
        )

ax_b.axhline(
    0,
    color="#555555",
    linewidth=1,
    linestyle="--"
)

ax_b.set_xticks(positions)
ax_b.set_xticklabels(
    labels,
    fontsize=12
)

ax_b.set_ylabel(
    "Change in mean programme score\n(24 h − 6 h)",
    fontsize=12
)

ax_b.set_title(
    "Temporal shifts in core Hallmark programmes",
    fontsize=15,
    fontweight="bold",
    pad=12
)

ax_b.grid(
    axis="y",
    color="#E5E5E5",
    linewidth=0.7,
    zorder=0
)

ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="-",
        color=BLUE,
        linewidth=3,
        markersize=7,
        label="Full"
    ),
    Line2D(
        [0], [0],
        marker="o",
        linestyle="-",
        color=ORANGE,
        linewidth=3,
        markersize=7,
        label="Landmark"
    )
]

ax_b.legend(
    handles=handles,
    frameon=False,
    ncol=2,
    loc="lower left",
    fontsize=11
)

fig.text(
    0.035, 0.54,
    "B",
    fontsize=20,
    fontweight="bold",
    va="top"
)


# ============================================================
# 7. METHODOLOGICAL NOTE — TWO LINES
# ============================================================

fig.text(
    0.5, 0.085,
    "Seven matched experimental blocks across A549, PC3 and VCAP. "
    "Points represent individual blocks; horizontal lines indicate means.",
    ha="center",
    fontsize=10.5
)

fig.text(
    0.5, 0.058,
    "Programme scores are mean Level 5 MODZ values, not GSEA NES.",
    ha="center",
    fontsize=10.5
)


# ============================================================
# 8. EXPORT
# ============================================================

pdf_path = OUTPUT / "Figure_6_temporal_sensitivity.pdf"
png_path = OUTPUT / "Figure_6_temporal_sensitivity.png"
csv_path = OUTPUT / "Figure_6_temporal_values.csv"

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

combined.to_csv(
    csv_path,
    index=False
)

plt.close(fig)

print("\nFIGURE 6 GENERATED SUCCESSFULLY")
print("PDF:", pdf_path)
print("PNG:", png_path)
print("CSV:", csv_path)
