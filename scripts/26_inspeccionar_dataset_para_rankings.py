import numpy as np
import pandas as pd
from pathlib import Path

print("=" * 80)
print("INSPECCIÓN DEL DATASET PRIMARIO PARA CONSTRUIR RANKINGS")
print("=" * 80)

# Buscar archivos NPZ
archivos_npz = list(Path(".").glob("*.npz"))

print("\nArchivos NPZ encontrados:")
for archivo in archivos_npz:
    print("-", archivo.name)

print("\n" + "=" * 80)

# Inspeccionar cada NPZ
for archivo in archivos_npz:
    print(f"\nARCHIVO: {archivo.name}")

    datos = np.load(archivo, allow_pickle=True)

    print("Claves:")
    for clave in datos.files:
        valor = datos[clave]
        print(f"  {clave}: shape = {valor.shape}, dtype = {valor.dtype}")

print("\n" + "=" * 80)

# Buscar archivos CSV relacionados con el dataset primario
archivos_csv = [
    p for p in Path(".").glob("*.csv")
    if "primario" in p.name.lower()
]

print("\nCSV relacionados con el dataset primario:")
for archivo in archivos_csv:
    print("-", archivo.name)

print("\nFIN DE LA INSPECCIÓN")