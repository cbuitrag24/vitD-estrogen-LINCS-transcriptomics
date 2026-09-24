import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

print("=" * 90)
print("FIGURA 4 — CONVERGENCIA FUNCIONAL HALLMARK")
print("=" * 90)

# ============================================================
# ARCHIVO PRINCIPAL
# ============================================================

archivo = Path("GSEA_FULL_vs_LANDMARK_comparacion.csv")

if not archivo.exists():
    raise FileNotFoundError(
        "No se encontró GSEA_FULL_vs_LANDMARK_comparacion.csv"
    )

df = pd.read_csv(archivo)

print("\nArchivo encontrado.")
print("Filas:", len(df))
print("Columnas:", list(df.columns))

# ============================================================
# HALLMARKS CENTRALES
# ============================================================

hallmarks = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1"
]

tratamientos = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

# Etiquetas para publicación
hallmark_labels = {
    "HALLMARK_E2F_TARGETS":
        "E2F TARGETS",

    "HALLMARK_G2M_CHECKPOINT":
        "G2M CHECKPOINT",

    "HALLMARK_MYC_TARGETS_V1":
        "MYC TARGETS V1"
}

treatment_labels = {
    "estradiol":
        "Estradiol",

    "calcitriol":
        "Calcitriol",

    "calcipotriol":
        "Calcipotriol"
}

# ============================================================
# FILTRAR DATOS
# ============================================================

core = df[
    df["hallmark"].isin(hallmarks)
    &
    df["tratamiento"].isin(tratamientos)
].copy()

print("\nDatos seleccionados para Figura 4:")
print(
    core[
        [
            "tratamiento",
            "hallmark",
            "NES_full",
            "q_full",
            "NES_landmark",
            "q_landmark",
            "robustez"
        ]
    ].sort_values(
        [
            "hallmark",
            "tratamiento"
        ]
    ).to_string(index=False)
)

# Deben ser exactamente 3 Hallmarks × 3 tratamientos
if len(core) != 9:
    raise ValueError(
        f"Se esperaban 9 filas y se encontraron {len(core)}."
    )

# ============================================================
# CONSTRUIR MATRICES
# ============================================================

nes_full = np.zeros((3, 3))
q_full = np.ones((3, 3))

nes_landmark = np.zeros((3, 3))
q_landmark = np.ones((3, 3))

for i, hallmark in enumerate(hallmarks):

    for j, tratamiento in enumerate(tratamientos):

        fila = core[
            (core["hallmark"] == hallmark)
            &
            (core["tratamiento"] == tratamiento)
        ]

        if len(fila) != 1:
            raise ValueError(
                f"Problema con {hallmark} / {tratamiento}"
            )

        nes_full[i, j] = float(
            fila["NES_full"].iloc[0]
        )

        q_full[i, j] = float(
            fila["q_full"].iloc[0]
        )

        nes_landmark[i, j] = float(
            fila["NES_landmark"].iloc[0]
        )

        q_landmark[i, j] = float(
            fila["q_landmark"].iloc[0]
        )

# ============================================================
# AUDITORÍA
# ============================================================

print("\nMatriz NES Full:")
print(nes_full)

print("\nMatriz NES Landmark:")
print(nes_landmark)

print("\nTodos los q Full < 0.05:",
      np.all(q_full < 0.05))

print("Todos los q Landmark < 0.05:",
      np.all(q_landmark < 0.05))

# ============================================================
# ESCALA COMÚN PARA LOS DOS HEATMAPS
# ============================================================

max_abs = max(
    np.abs(nes_full).max(),
    np.abs(nes_landmark).max()
)

limite = np.ceil(
    max_abs * 10
) / 10

norm = TwoSlopeNorm(
    vmin=-limite,
    vcenter=0,
    vmax=limite
)

# ============================================================
# FUNCIÓN PARA TEXTO DE SIGNIFICANCIA
# ============================================================

def estrellas(q):

    if q < 0.001:
        return "***"

    elif q < 0.01:
        return "**"

    elif q < 0.05:
        return "*"

    return ""

# ============================================================
# FIGURA
# ============================================================

fig, axes = plt.subplots(
    1,
    2,
    figsize=(10.8, 5.2),
    constrained_layout=True
)

# ============================================================
# PANEL A — FULL
# ============================================================

ax1 = axes[0]

im1 = ax1.imshow(
    nes_full,
    cmap="RdBu_r",
    norm=norm,
    aspect="auto"
)

ax1.set_title(
    "Full L1000 gene space",
    fontsize=11
)

ax1.set_xticks(
    np.arange(len(tratamientos))
)

ax1.set_xticklabels(
    [
        treatment_labels[x]
        for x in tratamientos
    ],
    fontsize=9
)

ax1.set_yticks(
    np.arange(len(hallmarks))
)

ax1.set_yticklabels(
    [
        hallmark_labels[x]
        for x in hallmarks
    ],
    fontsize=9
)

ax1.text(
    -0.16,
    1.07,
    "A",
    transform=ax1.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

# Texto dentro de las celdas
for i in range(3):

    for j in range(3):

        valor = nes_full[i, j]
        q = q_full[i, j]

        texto = (
    f"NES {valor:.2f}"
    f"{estrellas(q)}"
        )

        color_texto = (
            "white"
            if abs(valor) > limite * 0.48
            else "black"
        )

        ax1.text(
            j,
            i,
            texto,
            ha="center",
            va="center",
            fontsize=8.5,
            color=color_texto,
            fontweight="bold"
            if q < 0.05
            else "normal"
        )

# ============================================================
# PANEL B — LANDMARK
# ============================================================

ax2 = axes[1]

im2 = ax2.imshow(
    nes_landmark,
    cmap="RdBu_r",
    norm=norm,
    aspect="auto"
)

ax2.set_title(
    "Landmark genes",
    fontsize=11
)

ax2.set_xticks(
    np.arange(len(tratamientos))
)

ax2.set_xticklabels(
    [
        treatment_labels[x]
        for x in tratamientos
    ],
    fontsize=9
)

ax2.set_yticks(
    np.arange(len(hallmarks))
)

ax2.set_yticklabels(
    [
        hallmark_labels[x]
        for x in hallmarks
    ],
    fontsize=9
)

ax2.text(
    -0.16,
    1.07,
    "B",
    transform=ax2.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

for i in range(3):

    for j in range(3):

        valor = nes_landmark[i, j]
        q = q_landmark[i, j]

        texto = (
    f"NES {valor:.2f}"
    f"{estrellas(q)}"
        )

        color_texto = (
            "white"
            if abs(valor) > limite * 0.48
            else "black"
        )

        ax2.text(
            j,
            i,
            texto,
            ha="center",
            va="center",
            fontsize=8.5,
            color=color_texto,
            fontweight="bold"
            if q < 0.05
            else "normal"
        )

# ============================================================
# COLORBAR COMÚN
# ============================================================

cbar = fig.colorbar(
    im2,
    ax=axes,
    fraction=0.035,
    pad=0.04
)

cbar.set_label(
    "Normalised enrichment score (NES)",
    fontsize=9
)

# ============================================================
# NOTA ESTADÍSTICA
# ============================================================

fig.text(
    0.5,
    -0.025,
    "Negative NES indicates negative enrichment. "
    "* BH-adjusted empirical P < 0.05; ** P < 0.01; *** P < 0.001.",
    ha="center",
    va="top",
    fontsize=8.5
)

# ============================================================
# GUARDAR FIGURA
# ============================================================

fig.savefig(
    "Figure_4_Hallmark_core_convergence.pdf",
    bbox_inches="tight"
)

fig.savefig(
    "Figure_4_Hallmark_core_convergence.png",
    dpi=600,
    bbox_inches="tight"
)

# ============================================================
# GUARDAR DATOS EXACTOS UTILIZADOS
# ============================================================

tabla_figura4 = core[
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

tabla_figura4.to_csv(
    "Figure_4_Hallmark_core_values.csv",
    index=False
)

plt.close(fig)

# ============================================================
# FINAL
# ============================================================

print("\n" + "=" * 90)
print("ARCHIVOS GENERADOS")
print("=" * 90)

print(
    "- Figure_4_Hallmark_core_convergence.pdf"
)

print(
    "- Figure_4_Hallmark_core_convergence.png"
)

print(
    "- Figure_4_Hallmark_core_values.csv"
)

print("\nFIGURA 4 COMPLETADA")