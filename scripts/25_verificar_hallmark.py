from pathlib import Path

archivo = Path("data/h.all.v2026.1.Hs.symbols.gmt")

if not archivo.exists():
    print("ERROR: no se encontró el archivo Hallmark")
else:
    with open(archivo, "r", encoding="utf-8") as f:
        lineas = [linea.strip() for linea in f if linea.strip()]

    print("Archivo encontrado:", archivo.name)
    print("Número de gene sets Hallmark:", len(lineas))

    print("\nPrimeros 5 gene sets:")
    for linea in lineas[:5]:
        nombre = linea.split("\t")[0]
        print("-", nombre)