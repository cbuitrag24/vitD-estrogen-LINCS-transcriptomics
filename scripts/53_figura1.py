import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ============================================================
# FIGURA 1
# Study design and global transcriptional structure
# ============================================================

print("=" * 80)
print("FIGURA 1 — DISEÑO DEL ESTUDIO Y PCA")
print("=" * 80)

# ------------------------------------------------------------
# 1. Localizar archivos del análisis primario
# ------------------------------------------------------------

archivo_npz = Path("dataset_primario_expresion_12x12328.npz")

if not archivo_npz.exists():
    raise FileNotFoundError(
        "No se encontró dataset_primario_expresion_12x12328.npz"
    )

# Buscar automáticamente un archivo de metadatos primarios
candidatos_meta = [
    Path("dataset_primario_metadata.csv"),
    Path("metadata_primario.csv"),
    Path("diseno_primario.csv"),
    Path("firmas_primarias.csv"),
]

archivo_meta = None

for candidato in candidatos_meta:
    if candidato.exists():
        archivo_meta = candidato
        break

if archivo_meta is None:
    # Si el nombre es diferente, buscar CSV candidatos
    candidatos = list(Path(".").glob("*.csv"))

    posibles = [
        p for p in candidatos
        if (
            "primar" in p.name.lower()
            or "firma" in p.name.lower()
        )
    ]

    print("\nNo pude identificar automáticamente "
          "el archivo de metadatos primarios.")
    print("CSV candidatos encontrados:")

    for p in posibles:
        print(" -", p.name)

    raise FileNotFoundError(
        "Necesitamos identificar el CSV de metadatos "
        "correspondiente a las 12 firmas primarias."
    )

print("\nMatriz primaria:", archivo_npz)
print("Metadatos:", archivo_meta)

# ------------------------------------------------------------
# 2. Cargar matriz
# ------------------------------------------------------------

datos = np.load(archivo_npz)

if "expresion" not in datos.files:
    raise KeyError(
        "El NPZ no contiene la clave 'expresion'."
    )

X = datos["expresion"]

print("Dimensión de expresión:", X.shape)

if X.shape != (12, 12328):
    raise ValueError(
        f"Se esperaba una matriz (12, 12328), "
        f"pero se obtuvo {X.shape}."
    )

# ------------------------------------------------------------
# 3. Cargar metadatos
# ------------------------------------------------------------

meta = pd.read_csv(archivo_meta)

print("Filas de metadatos:", len(meta))
print("Columnas:", list(meta.columns))

if len(meta) != 12:
    raise ValueError(
        "El archivo de metadatos identificado no tiene 12 filas."
    )

# Identificar columnas
def buscar_columna(opciones):
    for opcion in opciones:
        for col in meta.columns:
            if col.lower() == opcion.lower():
                return col
    return None

col_cell = buscar_columna(
    ["cell_id", "cell", "cell_line", "celula"]
)

col_trat = buscar_columna(
    ["pert_iname", "treatment", "tratamiento", "compound"]
)

if col_cell is None or col_trat is None:
    raise ValueError(
        "No pude identificar las columnas de célula y tratamiento."
    )

print("\nColumna célula:", col_cell)
print("Columna tratamiento:", col_trat)

print("\nDiseño detectado:")
print(
    meta[[col_cell, col_trat]]
    .to_string(index=False)
)

# ------------------------------------------------------------
# 4. PCA
# ------------------------------------------------------------
#
# NO estandarizamos genes individualmente.
# Los valores MODZ ya representan respuestas diferenciales
# comparables y queremos reproducir la PCA realizada previamente.
# ------------------------------------------------------------

pca = PCA(n_components=2)

coords = pca.fit_transform(X)

varianza = pca.explained_variance_ratio_ * 100

print("\nVarianza explicada:")
print(f"PC1 = {varianza[0]:.2f}%")
print(f"PC2 = {varianza[1]:.2f}%")

pca_df = pd.DataFrame({
    "PC1": coords[:, 0],
    "PC2": coords[:, 1],
    "cell_id": meta[col_cell].astype(str),
    "treatment": meta[col_trat].astype(str)
})

# ------------------------------------------------------------
# 5. Normalizar nombres para presentación
# ------------------------------------------------------------

def nombre_tratamiento(x):
    s = str(x).lower()

    if "estradiol" in s:
        return "Estradiol"

    if "calcipotriol" in s:
        return "Calcipotriol"

    if "calcitriol" in s:
        return "Calcitriol"

    return str(x)

pca_df["treatment_label"] = (
    pca_df["treatment"]
    .apply(nombre_tratamiento)
)

# ------------------------------------------------------------
# 6. Crear Figura 1
# ------------------------------------------------------------

fig = plt.figure(
    figsize=(12, 5.5),
    constrained_layout=True
)

gs = fig.add_gridspec(
    1,
    2,
    width_ratios=[0.9, 1.25]
)

# ============================================================
# PANEL A — Diseño
# ============================================================

ax1 = fig.add_subplot(gs[0, 0])
ax1.axis("off")

ax1.text(
    0.02,
    0.97,
    "A",
    transform=ax1.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top"
)

ax1.text(
    0.50,
    0.92,
    "Matched LINCS L1000 study design",
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold"
)

# Caja 1: selección experimental
ax1.text(
    0.50,
    0.78,
   "LINCS L1000 perturbational signatures\n10 µM • 24 h",
    ha="center",
    va="center",
    fontsize=10,
    bbox=dict(
        boxstyle="round,pad=0.45",
        fill=False,
        linewidth=1
    )
)

# Flecha
ax1.annotate(
    "",
    xy=(0.50, 0.64),
    xytext=(0.50, 0.70),
    arrowprops=dict(
        arrowstyle="->",
        linewidth=1
    )
)

# Caja 2: contextos celulares
ax1.text(
    0.50,
    0.56,
    "4 matched experimental blocks\nA549 • HA1E • PC3 • VCAP",
    ha="center",
    va="center",
    fontsize=10,
    bbox=dict(
        boxstyle="round,pad=0.45",
        fill=False,
        linewidth=1
    )
)

# Flecha
ax1.annotate(
    "",
    xy=(0.50, 0.42),
    xytext=(0.50, 0.48),
    arrowprops=dict(
        arrowstyle="->",
        linewidth=1
    )
)

# Caja 3: tratamientos
ax1.text(
    0.50,
    0.34,
    "Estradiol  |  Calcitriol  |  Calcipotriol",
    ha="center",
    va="center",
    fontsize=10,
    bbox=dict(
        boxstyle="round,pad=0.45",
        fill=False,
        linewidth=1
    )
)

# Flecha
ax1.annotate(
    "",
    xy=(0.50, 0.20),
    xytext=(0.50, 0.26),
    arrowprops=dict(
        arrowstyle="->",
        linewidth=1
    )
)

# Caja 4: salida analítica
ax1.text(
    0.50,
    0.11,
    "12 Level 5 MODZ signatures\n4 blocks × 3 treatments",
    ha="center",
    va="center",
    fontsize=10,
    bbox=dict(
        boxstyle="round,pad=0.45",
        fill=False,
        linewidth=1
    )
)

ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)

# ============================================================
# PANEL B — PCA
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

trat_unicos = [
    "Estradiol",
    "Calcitriol",
    "Calcipotriol"
]

markers = {
    "Estradiol": "o",
    "Calcitriol": "s",
    "Calcipotriol": "^"
}

for tratamiento in trat_unicos:

    temp = pca_df[
        pca_df["treatment_label"] == tratamiento
    ]

    ax2.scatter(
        temp["PC1"],
        temp["PC2"],
        marker=markers[tratamiento],
        s=70,
        label=tratamiento,
        edgecolors="black",
        linewidths=0.6
    )

    for _, fila in temp.iterrows():

        ax2.annotate(
            fila["cell_id"],
            (
                fila["PC1"],
                fila["PC2"]
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8
        )

ax2.axhline(
    0,
    linewidth=0.5,
    alpha=0.5
)

ax2.axvline(
    0,
    linewidth=0.5,
    alpha=0.5
)

ax2.set_xlabel(
    f"PC1 ({varianza[0]:.1f}% variance)"
)

ax2.set_ylabel(
    f"PC2 ({varianza[1]:.1f}% variance)"
)

ax2.set_title(
    "Global transcriptional structure",
    fontsize=11
)

ax2.legend(
    frameon=False,
    fontsize=9
)

# ------------------------------------------------------------
# 7. Guardar
# ------------------------------------------------------------

fig.savefig(
    "Figure_1_study_design_PCA.pdf",
    bbox_inches="tight"
)

fig.savefig(
    "Figure_1_study_design_PCA.png",
    dpi=600,
    bbox_inches="tight"
)

pca_df.to_csv(
    "Figure_1_PCA_coordinates.csv",
    index=False
)

plt.close(fig)

print("\n" + "=" * 80)
print("ARCHIVOS GENERADOS")
print("=" * 80)

print("- Figure_1_study_design_PCA.pdf")
print("- Figure_1_study_design_PCA.png")
print("- Figure_1_PCA_coordinates.csv")

print("\nFIGURA 1 COMPLETADA")