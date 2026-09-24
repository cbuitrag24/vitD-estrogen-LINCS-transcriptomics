import numpy as np
import pandas as pd

# ============================================================
# VALIDACION FINAL DEL DATASET PRIMARIO
# LINCS L1000 - 24 h - 10 uM - diseño pareado
# ============================================================

ARCHIVO_EXPRESION = "dataset_primario_expresion_12x12328.npz"
ARCHIVO_METADATA = "dataset_primario_metadata.csv"
ARCHIVO_GENES = "dataset_primario_genes.csv"

print("=" * 70)
print("VALIDACION DEL DATASET PRIMARIO")
print("=" * 70)

# ------------------------------------------------------------
# 1. Cargar archivos
# ------------------------------------------------------------

npz = np.load(ARCHIVO_EXPRESION)
print("\nClaves NPZ:")
print(npz.files)

if "expresion" in npz.files:
    X = npz["expresion"]
elif "expression" in npz.files:
    X = npz["expression"]
else:
    raise ValueError("No se encontro la matriz de expresion en el NPZ.")

meta = pd.read_csv(ARCHIVO_METADATA)
genes = pd.read_csv(ARCHIVO_GENES)

print("\nForma matriz:", X.shape)
print("Firmas metadata:", len(meta))
print("Genes metadata:", len(genes))

# ------------------------------------------------------------
# 2. Controles de dimensiones
# ------------------------------------------------------------

if X.shape != (12, 12328):
    raise ValueError(
        f"Forma inesperada de la matriz: {X.shape}"
    )

if len(meta) != 12:
    raise ValueError(
        f"Se esperaban 12 firmas y hay {len(meta)}"
    )

if len(genes) != 12328:
    raise ValueError(
        f"Se esperaban 12328 genes y hay {len(genes)}"
    )

print("\nDimensiones: OK")

# ------------------------------------------------------------
# 3. Integridad numerica
# ------------------------------------------------------------

n_nan = int(np.isnan(X).sum())
n_inf = int(np.isinf(X).sum())

print("\nNaN:", n_nan)
print("Infinitos:", n_inf)

if n_nan != 0 or n_inf != 0:
    raise ValueError(
        "La matriz contiene valores no validos."
    )

print("Integridad numerica: OK")

# ------------------------------------------------------------
# 4. Comprobar tratamientos
# ------------------------------------------------------------

tratamientos_esperados = {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}

tratamientos = set(meta["pert_iname"].astype(str))

print("\nTratamientos encontrados:")
print(sorted(tratamientos))

if tratamientos != tratamientos_esperados:
    raise ValueError(
        "Los tratamientos no coinciden con los esperados."
    )

conteo_tratamientos = (
    meta["pert_iname"]
    .value_counts()
    .sort_index()
)

print("\nFirmas por tratamiento:")
print(conteo_tratamientos.to_string())

if not all(conteo_tratamientos == 4):
    raise ValueError(
        "No hay exactamente 4 firmas por tratamiento."
    )

print("Balance por tratamiento: OK")

# ------------------------------------------------------------
# 5. Comprobar lineas celulares
# ------------------------------------------------------------

celulas_esperadas = {
    "A549",
    "HA1E",
    "PC3",
    "VCAP",
}

celulas = set(meta["cell_id"].astype(str))

print("\nLineas celulares encontradas:")
print(sorted(celulas))

if celulas != celulas_esperadas:
    raise ValueError(
        "Las lineas celulares no coinciden con las esperadas."
    )

print("Lineas celulares: OK")

# ------------------------------------------------------------
# 6. Verificar diseño 1:1:1 por linea celular
# ------------------------------------------------------------

tabla = pd.crosstab(
    meta["cell_id"],
    meta["pert_iname"]
)

print("\nTabla linea celular x tratamiento:")
print(tabla.to_string())

if not (tabla == 1).all().all():
    raise ValueError(
        "El diseño no es 1:1:1 en todas las lineas celulares."
    )

print("\nDiseño 1:1:1: OK")

# ------------------------------------------------------------
# 7. Resumen numerico global
# ------------------------------------------------------------

print("\nResumen numerico de la matriz:")
print("Minimo:", float(np.min(X)))
print("Maximo:", float(np.max(X)))
print("Media:", float(np.mean(X)))
print("Mediana:", float(np.median(X)))
print("Desvio estandar:", float(np.std(X)))

# ------------------------------------------------------------
# 8. Genes landmark
# ------------------------------------------------------------

genes["pr_is_lm"] = pd.to_numeric(
    genes["pr_is_lm"],
    errors="coerce"
)

n_landmark = int((genes["pr_is_lm"] == 1).sum())

print("\nGenes landmark:", n_landmark)

if n_landmark != 978:
    raise ValueError(
        f"Se esperaban 978 genes landmark y hay {n_landmark}"
    )

print("Landmark genes: OK")

# ------------------------------------------------------------
# 9. Guardar resumen de validacion
# ------------------------------------------------------------

with open(
    "validacion_dataset_primario.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write("VALIDACION DEL DATASET PRIMARIO\n")
    f.write("=" * 50 + "\n\n")

    f.write(f"Forma matriz: {X.shape}\n")
    f.write(f"Firmas: {len(meta)}\n")
    f.write(f"Genes: {len(genes)}\n")
    f.write(f"Genes landmark: {n_landmark}\n")
    f.write(f"NaN: {n_nan}\n")
    f.write(f"Infinitos: {n_inf}\n")

    f.write("\nFirmas por tratamiento:\n")
    f.write(conteo_tratamientos.to_string())

    f.write("\n\nTabla linea celular x tratamiento:\n")
    f.write(tabla.to_string())

print("\nArchivo creado:")
print(" - validacion_dataset_primario.txt")

print("\nVALIDACION FINALIZADA CORRECTAMENTE")