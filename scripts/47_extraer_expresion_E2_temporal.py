import h5py
import numpy as np
import pandas as pd
from pathlib import Path

print("=" * 95)
print("EXTRACCIÓN DE EXPRESIÓN E2 TEMPORAL: 7 PARES, 6 H vs 24 H")
print("=" * 95)

# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_pares = Path("E2_temporal_firmas_pareadas.csv")

archivo_genes = Path(
    "data/GSE92742_Broad_LINCS_gene_info.txt.gz"
)

archivo_gctx = Path(
    r"D:\PAPER CB Y API NUEVO AGOSTO 2026"
    r"\GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
)

for archivo in [archivo_pares, archivo_genes, archivo_gctx]:
    if not archivo.exists():
        raise FileNotFoundError(
            f"No se encontró:\n{archivo}"
        )

print("\nArchivos encontrados correctamente.")

# ------------------------------------------------------------
# 2. Leer las 14 firmas seleccionadas
# ------------------------------------------------------------

meta = pd.read_csv(archivo_pares)

meta["tiempo"] = pd.to_numeric(
    meta["tiempo"],
    errors="raise"
).astype(int)

meta = meta.sort_values(
    ["cell_id", "bloque", "tiempo"]
).reset_index(drop=True)

sig_ids = meta["sig_id"].astype(str).tolist()

print("\nNúmero de firmas solicitadas:", len(sig_ids))

if len(sig_ids) != 14:
    raise ValueError(
        f"Esperábamos exactamente 14 firmas y encontramos {len(sig_ids)}."
    )

if meta["sig_id"].duplicated().any():
    raise ValueError("Hay sig_id duplicados en el archivo de pares.")

# Verificar estructura 7 pares x 2 tiempos
auditoria = (
    meta.groupby(["cell_id", "bloque"])["tiempo"]
    .apply(lambda x: sorted(x.tolist()))
)

if not all(x == [6, 24] for x in auditoria):
    raise ValueError(
        "Algún bloque no contiene exactamente los tiempos 6 h y 24 h."
    )

print("Estructura temporal verificada: 7 pares × 2 tiempos.")

# ------------------------------------------------------------
# 3. Leer información de genes
# ------------------------------------------------------------

genes = pd.read_csv(
    archivo_genes,
    sep="\t",
    low_memory=False
)

if len(genes) != 12328:
    raise ValueError(
        f"Esperábamos 12328 genes y encontramos {len(genes)}."
    )

print("Genes esperados:", len(genes))

# ------------------------------------------------------------
# 4. Abrir GCTX
# ------------------------------------------------------------

print("\nAbriendo GCTX...")
print("Esto puede tardar un poco. No cierres VS Code.")

with h5py.File(archivo_gctx, "r") as f:

    # Rutas estándar del GCTX LINCS
    matriz = f["0/DATA/0/matrix"]
    ids_columnas_raw = f["0/META/COL/id"][:]
    ids_filas_raw = f["0/META/ROW/id"][:]

    def decodificar(x):
        if isinstance(x, bytes):
            return x.decode("utf-8")
        return str(x)

    ids_columnas = [
        decodificar(x)
        for x in ids_columnas_raw
    ]

    ids_filas = [
        decodificar(x)
        for x in ids_filas_raw
    ]

    print("\nDimensión interna de la matriz:")
    print(matriz.shape)

    print("Número de IDs de columnas:", len(ids_columnas))
    print("Número de IDs de filas:", len(ids_filas))

    # --------------------------------------------------------
    # 5. Determinar orientación automáticamente
    # --------------------------------------------------------

    if (
        matriz.shape[0] == len(ids_filas)
        and matriz.shape[1] == len(ids_columnas)
    ):
        orientacion = "genes_x_firmas"

    elif (
        matriz.shape[0] == len(ids_columnas)
        and matriz.shape[1] == len(ids_filas)
    ):
        orientacion = "firmas_x_genes"

    else:
        raise ValueError(
            "No se pudo determinar la orientación de la matriz GCTX."
        )

    print("Orientación detectada:", orientacion)

    # --------------------------------------------------------
    # 6. Localizar las 14 firmas
    # --------------------------------------------------------

    mapa_columnas = {
        sig: i
        for i, sig in enumerate(ids_columnas)
    }

    faltantes = [
        sig
        for sig in sig_ids
        if sig not in mapa_columnas
    ]

    if faltantes:
        print("\nFIRMAS NO ENCONTRADAS:")
        for sig in faltantes:
            print(sig)

        raise ValueError(
            f"No se encontraron {len(faltantes)} firmas en el GCTX."
        )

    print("\nLas 14 firmas fueron encontradas en el GCTX.")

    indices = [
        mapa_columnas[sig]
        for sig in sig_ids
    ]

    # --------------------------------------------------------
    # 7. Extraer únicamente las 14 firmas
    # --------------------------------------------------------

    print("\nExtrayendo expresión...")

    if orientacion == "genes_x_firmas":

        # h5py exige índices ordenados
        orden = np.argsort(indices)
        indices_ordenados = np.array(indices)[orden]

        temp = matriz[:, indices_ordenados]

        # Recuperar el orden original solicitado
        inversa = np.argsort(orden)
        expresion = np.asarray(temp)[:, inversa].T

    else:

        orden = np.argsort(indices)
        indices_ordenados = np.array(indices)[orden]

        temp = matriz[indices_ordenados, :]

        inversa = np.argsort(orden)
        expresion = np.asarray(temp)[inversa, :]

# ------------------------------------------------------------
# 8. Comprobaciones
# ------------------------------------------------------------

print("\nDimensión extraída:", expresion.shape)

if expresion.shape != (14, 12328):
    raise ValueError(
        "La matriz extraída no tiene dimensión 14 × 12328."
    )

if not np.isfinite(expresion).all():
    raise ValueError(
        "La matriz contiene valores NaN o infinitos."
    )

print("Matriz verificada: 14 firmas × 12328 genes.")
print("Todos los valores son finitos.")

# ------------------------------------------------------------
# 9. Verificar orden de genes
# ------------------------------------------------------------

genes_gctx = [str(x) for x in ids_filas]

genes_info_ids = (
    genes["pr_gene_id"]
    .astype(str)
    .tolist()
)

if genes_gctx == genes_info_ids:
    print("Orden de genes verificado contra gene_info.")
else:
    print(
        "AVISO: el orden de genes del GCTX no coincide exactamente "
        "con gene_info. Se guardarán los IDs del GCTX para control."
    )

# ------------------------------------------------------------
# 10. Guardar dataset temporal
# ------------------------------------------------------------

np.savez_compressed(
    "dataset_E2_temporal_14x12328.npz",
    expresion=expresion,
    sig_ids=np.array(sig_ids),
    gene_ids=np.array(ids_filas)
)

meta.to_csv(
    "dataset_E2_temporal_metadata.csv",
    index=False
)

genes.to_csv(
    "dataset_E2_temporal_genes.csv",
    index=False
)

# ------------------------------------------------------------
# 11. Resumen
# ------------------------------------------------------------

print("\n" + "=" * 95)
print("RESUMEN FINAL")
print("=" * 95)

print("Firmas extraídas:", expresion.shape[0])
print("Genes:", expresion.shape[1])
print("Pares temporales:", len(meta) // 2)

print(
    "Células:",
    ", ".join(sorted(meta["cell_id"].unique()))
)

print("\nArchivos generados:")
print("- dataset_E2_temporal_14x12328.npz")
print("- dataset_E2_temporal_metadata.csv")
print("- dataset_E2_temporal_genes.csv")

print("\nEXTRACCIÓN COMPLETADA CORRECTAMENTE")