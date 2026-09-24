import numpy as np
import pandas as pd
from pathlib import Path

print("=" * 95)
print("CORRECCIÓN DEL ORDEN DE GENES — DATASET TEMPORAL E2")
print("=" * 95)

# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_npz = Path("dataset_E2_temporal_14x12328.npz")
archivo_gene_info = Path(
    "data/GSE92742_Broad_LINCS_gene_info.txt.gz"
)

if not archivo_npz.exists():
    raise FileNotFoundError(
        f"No se encontró: {archivo_npz}"
    )

if not archivo_gene_info.exists():
    raise FileNotFoundError(
        f"No se encontró: {archivo_gene_info}"
    )

# ------------------------------------------------------------
# 2. Leer dataset temporal
# ------------------------------------------------------------

datos = np.load(archivo_npz)

expresion = datos["expresion"]

gene_ids_gctx = [
    str(x) for x in datos["gene_ids"]
]

print("\nDimensión de expresión:", expresion.shape)
print("IDs de genes en GCTX:", len(gene_ids_gctx))

if expresion.shape != (14, 12328):
    raise ValueError(
        "La matriz temporal no tiene dimensión 14 × 12328."
    )

if len(gene_ids_gctx) != 12328:
    raise ValueError(
        "No hay exactamente 12328 IDs de genes."
    )

# ------------------------------------------------------------
# 3. Leer gene_info
# ------------------------------------------------------------

genes = pd.read_csv(
    archivo_gene_info,
    sep="\t",
    low_memory=False
)

genes["pr_gene_id"] = (
    genes["pr_gene_id"]
    .astype(str)
)

print("Filas en gene_info:", len(genes))

if len(genes) != 12328:
    raise ValueError(
        "gene_info no contiene exactamente 12328 genes."
    )

if genes["pr_gene_id"].duplicated().any():
    raise ValueError(
        "Hay IDs de genes duplicados en gene_info."
    )

# ------------------------------------------------------------
# 4. Verificar que sean los mismos genes
# ------------------------------------------------------------

set_gctx = set(gene_ids_gctx)
set_info = set(genes["pr_gene_id"])

if set_gctx != set_info:
    raise ValueError(
        "Los IDs del GCTX y gene_info no son idénticos."
    )

print("Los conjuntos de IDs son idénticos.")

# ------------------------------------------------------------
# 5. Reordenar gene_info según el GCTX
# ------------------------------------------------------------

genes_indexados = genes.set_index("pr_gene_id")

genes_ordenados = (
    genes_indexados
    .loc[gene_ids_gctx]
    .reset_index()
)

# ------------------------------------------------------------
# 6. Verificación estricta
# ------------------------------------------------------------

ids_finales = (
    genes_ordenados["pr_gene_id"]
    .astype(str)
    .tolist()
)

coincidencia = ids_finales == gene_ids_gctx

print(
    "\n¿Orden final idéntico al GCTX?",
    coincidencia
)

if not coincidencia:
    raise ValueError(
        "La corrección del orden no fue exitosa."
    )

if len(genes_ordenados) != expresion.shape[1]:
    raise ValueError(
        "El número de genes no coincide con "
        "las columnas de expresión."
    )

# ------------------------------------------------------------
# 7. Controles informativos
# ------------------------------------------------------------

print("\nPrimeros 10 genes después de corregir:")

print(
    genes_ordenados[
        [
            "pr_gene_id",
            "pr_gene_symbol",
            "pr_is_lm"
        ]
    ]
    .head(10)
    .to_string(index=False)
)

n_landmark = pd.to_numeric(
    genes_ordenados["pr_is_lm"],
    errors="coerce"
).fillna(0).astype(int).sum()

print("\nNúmero de genes Landmark:", n_landmark)

# ------------------------------------------------------------
# 8. Guardar tabla corregida
# ------------------------------------------------------------

salida = Path(
    "dataset_E2_temporal_genes_ordenados.csv"
)

genes_ordenados.to_csv(
    salida,
    index=False
)

# ------------------------------------------------------------
# 9. Resumen final
# ------------------------------------------------------------

print("\n" + "=" * 95)
print("RESUMEN FINAL")
print("=" * 95)

print("Firmas:", expresion.shape[0])
print("Genes:", expresion.shape[1])
print("Genes Landmark:", n_landmark)
print("Orden de genes verificado:", coincidencia)

print("\nArchivo generado:")
print("- dataset_E2_temporal_genes_ordenados.csv")

print(
    "\nIMPORTANTE: la matriz de expresión NO fue modificada."
)

print(
    "La tabla de genes fue reordenada para coincidir "
    "exactamente con las columnas de la matriz."
)

print("\nCORRECCIÓN COMPLETADA CORRECTAMENTE")