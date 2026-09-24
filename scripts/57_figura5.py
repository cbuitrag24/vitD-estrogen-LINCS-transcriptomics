import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

print("=" * 90)
print("FIGURA 5 — DEPENDENCIA DEL CONTEXTO CELULAR")
print("=" * 90)

# ============================================================
# ARCHIVOS
# ============================================================

archivo_full = Path("HALLMARK_CONTEXTO_matriz.csv")
archivo_lm = Path("HALLMARK_CONTEXTO_LM_consenso.csv")
archivo_resumen_full = Path("HALLMARK_CONTEXTO_resumen.csv")
archivo_resumen_lm = Path("HALLMARK_CONTEXTO_LM_resumen.csv")

for archivo in [
    archivo_full,
    archivo_lm,
    archivo_resumen_full,
    archivo_resumen_lm
]:
    if not archivo.exists():
        raise FileNotFoundError(
            f"No se encontró: {archivo.name}"
        )

full = pd.read_csv(archivo_full)
lm = pd.read_csv(archivo_lm)
resumen_full = pd.read_csv(archivo_resumen_full)
resumen_lm = pd.read_csv(archivo_resumen_lm)

print("\nArchivos cargados correctamente.")

# ============================================================
# ORDEN DE CONTEXTOS Y HALLMARKS
# ============================================================

contextos = [
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

hallmark_labels = {
    "HALLMARK_E2F_TARGETS":
        "E2F TARGETS",

    "HALLMARK_G2M_CHECKPOINT":
        "G2M CHECKPOINT",

    "HALLMARK_MYC_TARGETS_V1":
        "MYC TARGETS V1"
}

tratamientos = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

treatment_labels = {
    "estradiol": "Estradiol",
    "calcitriol": "Calcitriol",
    "calcipotriol": "Calcipotriol"
}

# ============================================================
# MATRIZ FULL
# ============================================================

filas_full = []

for cell in contextos:

    for hallmark in hallmarks:

        fila = full[
            (full["cell_id"] == cell)
            &
            (full["hallmark"] == hallmark)
        ]

        if len(fila) != 1:
            raise ValueError(
                f"Problema en Full con {cell} / {hallmark}"
            )

        for tratamiento in tratamientos:

            filas_full.append({
                "cell_id": cell,
                "hallmark": hallmark,
                "tratamiento": tratamiento,
                "score": float(
                    fila[tratamiento].iloc[0]
                )
            })

full_long = pd.DataFrame(filas_full)

# ============================================================
# MATRIZ LANDMARK
# ============================================================

mapa_lm_cols = {
    "estradiol": "estradiol_LM",
    "calcitriol": "calcitriol_LM",
    "calcipotriol": "calcipotriol_LM"
}

filas_lm = []

for cell in contextos:

    for hallmark in hallmarks:

        fila = lm[
            (lm["cell_id"] == cell)
            &
            (lm["hallmark"] == hallmark)
        ]

        if len(fila) != 1:
            raise ValueError(
                f"Problema en Landmark con {cell} / {hallmark}"
            )

        for tratamiento in tratamientos:

            col = mapa_lm_cols[tratamiento]

            filas_lm.append({
                "cell_id": cell,
                "hallmark": hallmark,
                "tratamiento": tratamiento,
                "score": float(
                    fila[col].iloc[0]
                )
            })

lm_long = pd.DataFrame(filas_lm)

# ============================================================
# AUDITORÍA DE RESÚMENES
# ============================================================

print("\nResumen Full:")
print(
    resumen_full.to_string(index=False)
)

print("\nResumen Landmark:")
print(
    resumen_lm.to_string(index=False)
)

# ============================================================
# MATRICES PARA HEATMAP
# ============================================================

# Filas = contexto × Hallmark
row_labels = []

for cell in contextos:
    for hallmark in hallmarks:
        row_labels.append(
            f"{cell} | {hallmark_labels[hallmark]}"
        )

def construir_matriz(df_long):

    matriz = np.zeros(
        (
            len(row_labels),
            len(tratamientos)
        )
    )

    fila_idx = 0

    for cell in contextos:

        for hallmark in hallmarks:

            subset = df_long[
                (df_long["cell_id"] == cell)
                &
                (df_long["hallmark"] == hallmark)
            ]

            for j, tratamiento in enumerate(tratamientos):

                fila = subset[
                    subset["tratamiento"]
                    == tratamiento
                ]

                if len(fila) != 1:
                    raise ValueError(
                        f"Problema en matriz: "
                        f"{cell}, {hallmark}, {tratamiento}"
                    )

                matriz[
                    fila_idx,
                    j
                ] = float(
                    fila["score"].iloc[0]
                )

            fila_idx += 1

    return matriz

matriz_full = construir_matriz(full_long)
matriz_lm = construir_matriz(lm_long)

print("\nMatriz Full:")
print(matriz_full)

print("\nMatriz Landmark:")
print(matriz_lm)

# ============================================================
# ESCALA COMÚN
# ============================================================

max_abs = max(
    np.abs(matriz_full).max(),
    np.abs(matriz_lm).max()
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
# FIGURA
# ============================================================

fig = plt.figure(
    figsize=(14.0, 10.5),
    constrained_layout=True
)

gs = fig.add_gridspec(
    2,
    2,
    height_ratios=[1.6, 0.8]
)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])

# ============================================================
# PANEL A — FULL
# ============================================================

im1 = ax1.imshow(
    matriz_full,
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
        treatment_labels[t]
        for t in tratamientos
    ],
    fontsize=9
)

ax1.set_yticks(
    np.arange(len(row_labels))
)

ax1.set_yticklabels(
    row_labels,
    fontsize=8
)

ax1.text(
    -0.20,
    1.03,
    "A",
    transform=ax1.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

for i in range(
    matriz_full.shape[0]
):

    for j in range(
        matriz_full.shape[1]
    ):

        valor = matriz_full[i, j]

        color_texto = (
            "white"
            if abs(valor) > limite * 0.45
            else "black"
        )

        ax1.text(
            j,
            i,
            f"{valor:.2f}",
            ha="center",
            va="center",
            fontsize=7.5,
            color=color_texto
        )

# ============================================================
# PANEL B — LANDMARK
# ============================================================

im2 = ax2.imshow(
    matriz_lm,
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
        treatment_labels[t]
        for t in tratamientos
    ],
    fontsize=9
)

ax2.set_yticks(
    np.arange(len(row_labels))
)

ax2.set_yticklabels(
    row_labels,
    fontsize=8
)

ax2.text(
    -0.20,
    1.03,
    "B",
    transform=ax2.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

for i in range(
    matriz_lm.shape[0]
):

    for j in range(
        matriz_lm.shape[1]
    ):

        valor = matriz_lm[i, j]

        color_texto = (
            "white"
            if abs(valor) > limite * 0.45
            else "black"
        )

        ax2.text(
            j,
            i,
            f"{valor:.2f}",
            ha="center",
            va="center",
            fontsize=7.5,
            color=color_texto
        )

# ============================================================
# COLORBAR COMÚN
# ============================================================

cbar = fig.colorbar(
    im2,
    ax=[ax1, ax2],
    fraction=0.028,
    pad=0.03
)

cbar.set_label(
    "Mean Hallmark programme score",
    fontsize=9
)

# ============================================================
# PANEL C — RESUMEN DE CONVERGENCIA POR CONTEXTO
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

res_full = (
    resumen_full
    .set_index("cell_id")
    .loc[contextos]
)

res_lm = (
    resumen_lm
    .set_index("cell_id")
    .loc[contextos]
)

x = np.arange(
    len(contextos)
)

width = 0.34

val_full = (
    res_full[
        "fraccion_hallmarks_negativos_3de3"
    ].values
)

val_lm = (
    res_lm[
        "fraccion_hallmarks_negativos_3de3"
    ].values
)

bars_full = ax3.bar(
    x - width / 2,
    val_full,
    width,
    label="Full"
)

bars_lm = ax3.bar(
    x + width / 2,
    val_lm,
    width,
    label="Landmark"
)

ax3.set_xticks(x)

ax3.set_xticklabels(
    contextos,
    fontsize=9
)

ax3.set_ylim(
    0,
    1.12
)

ax3.set_ylabel(
    "Fraction of core Hallmarks\nnegative in all 3 treatments"
)

ax3.set_title(
    "Context-level convergence of the core programmes",
    fontsize=11
)

ax3.legend(
    frameon=False,
    fontsize=9
)

# Etiquetas 0/3, 1/3, 3/3
for bars, valores in [
    (bars_full, val_full),
    (bars_lm, val_lm)
]:

    for bar, valor in zip(
        bars,
        valores
    ):

        n = int(
            round(
                valor * 3
            )
        )

        ax3.text(
            bar.get_x()
            + bar.get_width() / 2,
            bar.get_height()
            + 0.035,
            f"{n}/3",
            ha="center",
            va="bottom",
            fontsize=8.5
        )

ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

# ============================================================
# NOTA
# ============================================================

fig.text(
    0.5,
    -0.025,
    "Context-level values are descriptive because cell line and experimental block are confounded "
    "in the primary matched design.",
    ha="center",
    va="top",
    fontsize=8.5
)

# ============================================================
# GUARDAR
# ============================================================

fig.savefig(
    "Figure_5_context_dependence.pdf",
    bbox_inches="tight"
)

fig.savefig(
    "Figure_5_context_dependence.png",
    dpi=600,
    bbox_inches="tight"
)

# Datos utilizados
tabla_figura5 = pd.concat(
    [
        full_long.assign(
            analysis_space="Full"
        ),
        lm_long.assign(
            analysis_space="Landmark"
        )
    ],
    ignore_index=True
)

tabla_figura5.to_csv(
    "Figure_5_context_values.csv",
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
    "- Figure_5_context_dependence.pdf"
)

print(
    "- Figure_5_context_dependence.png"
)

print(
    "- Figure_5_context_values.csv"
)

print("\nFIGURA 5 COMPLETADA")