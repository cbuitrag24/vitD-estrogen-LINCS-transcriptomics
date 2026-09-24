import h5py
import numpy as np
import pandas as pd
from pathlib import Path

# ============================================================
# EXTRACCION DE FIRMAS LEVEL 5 - LINCS L1000
# ============================================================

# Archivo GCTX original
GCTX = Path("data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")

# Firmas seleccionadas previamente
ARCHIVO_FIRMAS = Path("data/curation/firmas_para_extraccion_level5.csv")

# Información de genes
ARCHIVO_GENES = Path(
    "data/GSE92742_Broad_LINCS_gene_info.txt.gz"
)

print("=" * 60)
print("EXTRACCION LEVEL 5")
print("=" * 60)

# ------------------------------------------------------------
# 1. Leer las firmas seleccionadas
# ------------------------------------------------------------

firmas = pd.read_csv(ARCHIVO_FIRMAS)

if "sig_id" not in firmas.columns:
    raise ValueError("No se encontro la columna sig_id.")

sig_seleccionados = firmas["sig_id"].astype(str).tolist()

print("\nFirmas solicitadas:", len(sig_seleccionados))
print("sig_id unicos:", len(set(sig_seleccionados)))

if len(sig_seleccionados) != len(set(sig_seleccionados)):
    raise ValueError("Hay sig_id duplicados en el archivo de seleccion.")

# ------------------------------------------------------------
# 2. Abrir GCTX
# ------------------------------------------------------------

print("\nAbriendo GCTX...")
print(GCTX)

with h5py.File(GCTX, "r") as f:

    matrix = f["0"]["DATA"]["0"]["matrix"]
    col_ids_raw = f["0"]["META"]["COL"]["id"]
    row_ids_raw = f["0"]["META"]["ROW"]["id"]

    print("\nForma matriz:", matrix.shape)
    print("Tipo:", matrix.dtype)

    # COL = firmas
    col_ids = [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in col_ids_raw[:]
    ]

    # ROW = genes
    row_ids = [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in row_ids_raw[:]
    ]

    print("Firmas en GCTX:", len(col_ids))
    print("Genes en GCTX:", len(row_ids))

    # --------------------------------------------------------
    # 3. Localizar los sig_id seleccionados
    # --------------------------------------------------------

    mapa_sig = {sig: i for i, sig in enumerate(col_ids)}

    faltantes = [
        sig for sig in sig_seleccionados
        if sig not in mapa_sig
    ]

    print("\nFirmas encontradas:",
          len(sig_seleccionados) - len(faltantes))

    print("Firmas faltantes:", len(faltantes))

    if faltantes:
        pd.DataFrame(
            {"sig_id_faltante": faltantes}
        ).to_csv(
            "firmas_faltantes_level5.csv",
            index=False
        )

        raise ValueError(
            "Hay firmas que no se encontraron en el GCTX. "
            "Se creo firmas_faltantes_level5.csv"
        )

    # Índices correspondientes a nuestras firmas
    indices_originales = np.array(
        [mapa_sig[sig] for sig in sig_seleccionados],
        dtype=int
    )

    # h5py requiere índices crecientes para lectura eficiente
    orden = np.argsort(indices_originales)
    indices_ordenados = indices_originales[orden]

    print("\nExtrayendo valores de expresion...")

    # La matriz está almacenada como firmas x genes
    expresion_ordenada = matrix[indices_ordenados, :]

    # Restaurar el orden original de firmas_para_extraccion_level5.csv
    inverso = np.argsort(orden)
    expresion = expresion_ordenada[inverso, :]

# ------------------------------------------------------------
# 4. Controles de dimensiones
# ------------------------------------------------------------

print("\nForma de la matriz extraida:", expresion.shape)

if expresion.shape[0] != len(sig_seleccionados):
    raise ValueError(
        "El numero de firmas extraidas no coincide."
    )

if expresion.shape[1] != len(row_ids):
    raise ValueError(
        "El numero de genes extraidos no coincide."
    )

# ------------------------------------------------------------
# 5. Guardar matriz numérica
# ------------------------------------------------------------

np.savez_compressed(
    "expresion_level5_213x12328.npz",
    expression=expresion,
    sig_ids=np.array(sig_seleccionados),
    gene_ids=np.array(row_ids)
)

# ------------------------------------------------------------
# 6. Guardar metadatos de las firmas
# ------------------------------------------------------------

firmas.to_csv(
    "metadata_firmas_level5.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Vincular IDs de genes con símbolos
# ------------------------------------------------------------

gene_info = pd.read_csv(
    ARCHIVO_GENES,
    sep="\t",
    dtype={"pr_gene_id": str}
)

genes_gctx = pd.DataFrame(
    {
        "matrix_position": np.arange(len(row_ids)),
        "pr_gene_id": row_ids
    }
)

genes_anotados = genes_gctx.merge(
    gene_info,
    on="pr_gene_id",
    how="left"
)

genes_anotados.to_csv(
    "metadata_genes_level5.csv",
    index=False
)

# ------------------------------------------------------------
# 8. Control de valores
# ------------------------------------------------------------

print("\nValores no finitos:",
      np.size(expresion) - np.isfinite(expresion).sum())

print("Valor minimo:", float(np.nanmin(expresion)))
print("Valor maximo:", float(np.nanmax(expresion)))
print("Media:", float(np.nanmean(expresion)))

print("\nGenes con simbolo:",
      genes_anotados["pr_gene_symbol"].notna().sum())

print("Genes sin simbolo:",
      genes_anotados["pr_gene_symbol"].isna().sum())

print("\nArchivos creados:")
print(" - expresion_level5_213x12328.npz")
print(" - metadata_firmas_level5.csv")
print(" - metadata_genes_level5.csv")

print("\nEXTRACCION FINALIZADA CORRECTAMENTE")