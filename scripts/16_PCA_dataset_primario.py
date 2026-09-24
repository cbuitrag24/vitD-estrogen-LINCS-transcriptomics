import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ============================================================
# PCA EXPLORATORIA DEL DATASET PRIMARIO
# LINCS L1000 - 24 h - 10 uM - diseño pareado
#
# Objetivo:
# explorar la estructura global de las 12 firmas
# sin realizar inferencia estadistica.
# ============================================================

ARCHIVO_EXPRESION = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_METADATA = "dataset_primario_metadata.csv"
ARCHIVO_GENES = "dataset_primario_genes.csv"

print("=" * 70)
print("PCA EXPLORATORIA - DATASET PRIMARIO")
print("=" * 70)

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------

npz = np.load(ARCHIVO_EXPRESION)

if "expresion" in npz.files:
    X = npz["expresion"]
elif "expression" in npz.files:
    X = npz["expression"]
else:
    raise ValueError("No se encontro la matriz de expresion.")

meta = pd.read_csv(ARCHIVO_METADATA)
genes = pd.read_csv(ARCHIVO_GENES)

print("\nForma matriz:")
print(X.shape)

# ------------------------------------------------------------
# 2. Control previo
# ------------------------------------------------------------

if X.shape[0] != len(meta):
    raise ValueError(
        "Las filas de expresion no coinciden con metadata."
    )

if X.shape[1] != len(genes):
    raise ValueError(
        "Las columnas de expresion no coinciden con genes."
    )

# ------------------------------------------------------------
# 3. PCA sobre los valores Level 5
# ------------------------------------------------------------
# Los perfiles Level 5 ya representan valores normalizados
# de perturbacion. Para esta PCA primaria NO escalamos
# individualmente cada gen, porque hacerlo daria el mismo peso
# a genes con variaciones biologicas muy diferentes.

pca = PCA(n_components=min(X.shape[0], X.shape[1]))

coords = pca.fit_transform(X)

varianza = pca.explained_variance_ratio_ * 100

print("\nVarianza explicada:")
for i in range(min(5, len(varianza))):
    print(
        f"PC{i+1}: {varianza[i]:.2f}%"
    )

# ------------------------------------------------------------
# 4. Guardar coordenadas
# ------------------------------------------------------------

resultado = meta.copy()

for i in range(coords.shape[1]):
    resultado[f"PC{i+1}"] = coords[:, i]

resultado.to_csv(
    "PCA_dataset_primario_coordenadas.csv",
    index=False
)

# ------------------------------------------------------------
# 5. Guardar varianza explicada
# ------------------------------------------------------------

tabla_varianza = pd.DataFrame(
    {
        "PC": [
            f"PC{i+1}"
            for i in range(len(varianza))
        ],
        "variance_percent": varianza,
        "cumulative_percent": np.cumsum(varianza),
    }
)

tabla_varianza.to_csv(
    "PCA_dataset_primario_varianza.csv",
    index=False
)

# ------------------------------------------------------------
# 6. Figura PCA: tratamiento
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 6))

tratamientos = [
    "estradiol",
    "calcitriol",
    "calcipotriol",
]

marcadores = {
    "estradiol": "o",
    "calcitriol": "s",
    "calcipotriol": "^",
}

for tratamiento in tratamientos:

    mask = (
        resultado["pert_iname"] == tratamiento
    )

    ax.scatter(
        resultado.loc[mask, "PC1"],
        resultado.loc[mask, "PC2"],
        marker=marcadores[tratamiento],
        s=90,
        label=tratamiento,
        alpha=0.8,
    )

# Etiquetar cada punto con linea celular
for _, fila in resultado.iterrows():

    ax.annotate(
        fila["cell_id"],
        (fila["PC1"], fila["PC2"]),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=8,
    )

ax.axhline(0, linewidth=0.5)
ax.axvline(0, linewidth=0.5)

ax.set_xlabel(
    f"PC1 ({varianza[0]:.1f}% variance)"
)

ax.set_ylabel(
    f"PC2 ({varianza[1]:.1f}% variance)"
)

ax.set_title(
    "PCA of paired LINCS L1000 signatures\n"
    "24 h, 10 µM"
)

ax.legend(
    title="Perturbagen",
    frameon=False
)

plt.tight_layout()

plt.savefig(
    "PCA_dataset_primario_tratamientos.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# 7. PCA usando SOLO genes landmark
# ------------------------------------------------------------

genes["pr_is_lm"] = pd.to_numeric(
    genes["pr_is_lm"],
    errors="coerce"
)

mask_landmark = (
    genes["pr_is_lm"].to_numpy() == 1
)

X_landmark = X[:, mask_landmark]

print("\nForma matriz landmark:")
print(X_landmark.shape)

pca_lm = PCA(
    n_components=min(
        X_landmark.shape[0],
        X_landmark.shape[1]
    )
)

coords_lm = pca_lm.fit_transform(
    X_landmark
)

varianza_lm = (
    pca_lm.explained_variance_ratio_
    * 100
)

resultado_lm = meta.copy()

resultado_lm["PC1"] = coords_lm[:, 0]
resultado_lm["PC2"] = coords_lm[:, 1]

resultado_lm.to_csv(
    "PCA_landmark_coordenadas.csv",
    index=False
)

tabla_varianza_lm = pd.DataFrame(
    {
        "PC": [
            f"PC{i+1}"
            for i in range(len(varianza_lm))
        ],
        "variance_percent": varianza_lm,
        "cumulative_percent": np.cumsum(
            varianza_lm
        ),
    }
)

tabla_varianza_lm.to_csv(
    "PCA_landmark_varianza.csv",
    index=False
)

# ------------------------------------------------------------
# 8. Figura landmark
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 6))

for tratamiento in tratamientos:

    mask = (
        resultado_lm["pert_iname"]
        == tratamiento
    )

    ax.scatter(
        resultado_lm.loc[mask, "PC1"],
        resultado_lm.loc[mask, "PC2"],
        marker=marcadores[tratamiento],
        s=90,
        label=tratamiento,
        alpha=0.8,
    )

for _, fila in resultado_lm.iterrows():

    ax.annotate(
        fila["cell_id"],
        (fila["PC1"], fila["PC2"]),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=8,
    )

ax.axhline(0, linewidth=0.5)
ax.axvline(0, linewidth=0.5)

ax.set_xlabel(
    f"PC1 ({varianza_lm[0]:.1f}% variance)"
)

ax.set_ylabel(
    f"PC2 ({varianza_lm[1]:.1f}% variance)"
)

ax.set_title(
    "PCA of landmark-gene LINCS signatures\n"
    "24 h, 10 µM"
)

ax.legend(
    title="Perturbagen",
    frameon=False
)

plt.tight_layout()

plt.savefig(
    "PCA_dataset_primario_landmark.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# 9. Resumen
# ------------------------------------------------------------

print("\nArchivos creados:")
print(" - PCA_dataset_primario_coordenadas.csv")
print(" - PCA_dataset_primario_varianza.csv")
print(" - PCA_dataset_primario_tratamientos.png")
print(" - PCA_landmark_coordenadas.csv")
print(" - PCA_landmark_varianza.csv")
print(" - PCA_dataset_primario_landmark.png")

print("\nPCA FINALIZADA CORRECTAMENTE")