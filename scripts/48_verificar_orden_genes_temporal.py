import numpy as np
import pandas as pd
from pathlib import Path

print("=" * 95)
print("VERIFICACIÓN DEL ORDEN DE GENES — DATASET TEMPORAL E2")
print("=" * 95)

# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_npz = Path("dataset_E2_temporal_14x12328.npz")
archivo_gene_info = Path("data/GSE92742_Broad_LINCS_gene_info.txt.gz")

if not archivo_npz.exists():
    raise FileNotFoundError(f"No se encontró: {archivo_npz}")

if not archivo_gene_info.exists():
    raise FileNotFoundError(f"No se encontró: {archivo_gene_info}")

# ------------------------------------------------------------
# 2. Leer IDs guardados directamente desde el GCTX
# ------------------------------------------------------------

datos = np.load(archivo_npz)

gene_ids_gctx = [
    str(x) for x in datos["gene_ids"]
]

print("\nGenes guardados desde GCTX:", len(gene_ids_gctx))

# ------------------------------------------------------------
# 3. Leer gene_info
# ------------------------------------------------------------

genes = pd.read_csv(
    archivo_gene_info,
    sep="\t",
    low_memory=False
)

gene_ids_info = (
    genes["pr_gene_id"]
    .astype(str)
    .tolist()
)

print("Genes en gene_info:", len(gene_ids_info))

# ------------------------------------------------------------
# 4. Comparación directa
# ------------------------------------------------------------

coinciden_directamente = (
    gene_ids_gctx == gene_ids_info
)

print("\n¿Coinciden exactamente y en el mismo orden?")
print(coinciden_directamente)

# ------------------------------------------------------------
# 5. Comparar como conjuntos
# ------------------------------------------------------------

set_gctx = set(gene_ids_gctx)
set_info = set(gene_ids_info)

solo_gctx = sorted(set_gctx - set_info)
solo_info = sorted(set_info - set_gctx)

print("\nIDs presentes sólo en GCTX:", len(solo_gctx))
print("IDs presentes sólo en gene_info:", len(solo_info))

# ------------------------------------------------------------
# 6. Buscar primera discrepancia
# ------------------------------------------------------------

primera = None

for i, (a, b) in enumerate(
    zip(gene_ids_gctx, gene_ids_info)
):
    if a != b:
        primera = i
        break

if primera is None:
    print("\nNo se encontraron discrepancias posición por posición.")
else:
    print("\nPrimera discrepancia en posición:", primera)
    print("ID GCTX     :", repr(gene_ids_gctx[primera]))
    print("ID gene_info:", repr(gene_ids_info[primera]))

# ------------------------------------------------------------
# 7. Mostrar primeros IDs para diagnóstico
# ------------------------------------------------------------

print("\nPrimeros 10 IDs del GCTX:")
for x in gene_ids_gctx[:10]:
    print(repr(x))

print("\nPrimeros 10 IDs de gene_info:")
for x in gene_ids_info[:10]:
    print(repr(x))

# ------------------------------------------------------------
# 8. Comprobar si el problema es sólo el orden
# ------------------------------------------------------------

if set_gctx == set_info:

    print("\n" + "=" * 95)
    print("RESULTADO")
    print("=" * 95)

    if coinciden_directamente:
        print("Los IDs y el orden coinciden exactamente.")
    else:
        print("Los mismos 12328 IDs están presentes en ambos archivos.")
        print("La diferencia corresponde únicamente al ORDEN.")
        print("Esto puede corregirse de forma segura mediante los IDs.")

else:

    print("\n" + "=" * 95)
    print("RESULTADO")
    print("=" * 95)

    print("Los conjuntos de IDs NO son idénticos.")

    if solo_gctx:
        print("\nEjemplos presentes sólo en GCTX:")
        print(solo_gctx[:10])

    if solo_info:
        print("\nEjemplos presentes sólo en gene_info:")
        print(solo_info[:10])

    print(
        "\nNo se debe continuar con el análisis hasta "
        "resolver esta diferencia."
    )

print("\nFIN")