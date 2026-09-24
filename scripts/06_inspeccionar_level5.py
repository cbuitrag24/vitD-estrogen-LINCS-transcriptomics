import gzip
import shutil
from pathlib import Path

archivo_gz = Path(
    "data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
)

archivo_gctx = archivo_gz.with_suffix("")

print("Archivo comprimido:")
print(archivo_gz)

print("\nExiste:", archivo_gz.exists())

if archivo_gz.exists():
    tamano_gb = archivo_gz.stat().st_size / (1024**3)
    print(f"Tamaño: {tamano_gb:.2f} GB")

print("\nTodavía NO vamos a descomprimirlo.")
print("Comprobación terminada.")