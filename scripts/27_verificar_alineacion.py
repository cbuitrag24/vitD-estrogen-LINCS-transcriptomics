import numpy as np
import pandas as pd

print("=" * 80)
print("VERIFICACIÓN DE ALINEACIÓN DEL DATASET PRIMARIO")
print("=" * 80)

# Cargar matriz de expresión
datos = np.load(
    "dataset_primario_expresion_12x12328.npz",
    allow_pickle=True
)
expresion = datos["expresion"]

# Cargar metadatos y genes
metadata = pd.read_csv("dataset_primario_metadata.csv")
genes = pd.read_csv("dataset_primario_genes.csv")

print("\nDimensiones:")
print("Matriz de expresión:", expresion.shape)
print("Metadata:", metadata.shape)
print("Genes:", genes.shape)

print("\nColumnas de metadata:")
print(metadata.columns.tolist())

print("\nColumnas de genes:")
print(genes.columns.tolist())

print("\nPrimeras 12 filas de metadata:")
print(metadata.to_string(index=False))

print("\nPrimeros 5 genes:")
print(genes.head().to_string(index=False))

print("\nControles:")
print("Filas expresión == filas metadata:",
      expresion.shape[0] == len(metadata))

print("Columnas expresión == filas genes:",
      expresion.shape[1] == len(genes))

print("\nFIN DE LA VERIFICACIÓN")