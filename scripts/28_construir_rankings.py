import numpy as np
import pandas as pd

print("=" * 80)
print("CONSTRUCCIÓN DE RANKINGS TRANSCRIPCIONALES")
print("=" * 80)

# ---------------------------------------------------------------------
# 1. Cargar datos
# ---------------------------------------------------------------------

datos = np.load(
    "dataset_primario_expresion_12x12328.npz",
    allow_pickle=True
)

expresion = datos["expresion"]

metadata = pd.read_csv("dataset_primario_metadata.csv")
genes = pd.read_csv("dataset_primario_genes.csv")

# ---------------------------------------------------------------------
# 2. Controles básicos
# ---------------------------------------------------------------------

if expresion.shape[0] != len(metadata):
    raise ValueError("ERROR: las filas de expresión no coinciden con metadata")

if expresion.shape[1] != len(genes):
    raise ValueError("ERROR: las columnas de expresión no coinciden con genes")

if genes["pr_gene_symbol"].isna().any():
    raise ValueError("ERROR: existen genes sin símbolo génico")

print("\nDimensión de la matriz:", expresion.shape)

# ---------------------------------------------------------------------
# 3. Tratamientos primarios
# ---------------------------------------------------------------------

tratamientos = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

print("\nNúmero de firmas por tratamiento:")

for tratamiento in tratamientos:

    mascara = (
        metadata["pert_iname"]
        .astype(str)
        .str.lower()
        .eq(tratamiento)
    )

    n = int(mascara.sum())

    print(f"- {tratamiento}: {n}")

    if n != 4:
        raise ValueError(
            f"ERROR: {tratamiento} debería tener exactamente 4 firmas, "
            f"pero se encontraron {n}"
        )

# ---------------------------------------------------------------------
# 4. Construir rankings
# ---------------------------------------------------------------------

resumen = []

for tratamiento in tratamientos:

    mascara = (
        metadata["pert_iname"]
        .astype(str)
        .str.lower()
        .eq(tratamiento)
    )

    indices = np.where(mascara.to_numpy())[0]

    matriz_tratamiento = expresion[indices, :]

    # Promedio del MODZ entre los cuatro contextos pareados
    score_medio = np.mean(matriz_tratamiento, axis=0)

    # Número de contextos con respuesta positiva y negativa
    n_positivo = np.sum(matriz_tratamiento > 0, axis=0)
    n_negativo = np.sum(matriz_tratamiento < 0, axis=0)

    ranking = pd.DataFrame({
        "pr_gene_id": genes["pr_gene_id"],
        "pr_gene_symbol": genes["pr_gene_symbol"],
        "pr_is_lm": genes["pr_is_lm"],
        "pr_is_bing": genes["pr_is_bing"],
        "mean_modz": score_medio,
        "n_contextos_positivo": n_positivo,
        "n_contextos_negativo": n_negativo
    })

    ranking = ranking.sort_values(
        "mean_modz",
        ascending=False
    ).reset_index(drop=True)

    archivo_salida = f"RANKING_{tratamiento}.csv"

    ranking.to_csv(
        archivo_salida,
        index=False
    )

    resumen.append({
        "tratamiento": tratamiento,
        "n_firmas": len(indices),
        "n_genes": len(ranking),
        "mean_modz_max": ranking["mean_modz"].max(),
        "mean_modz_min": ranking["mean_modz"].min()
    })

    print(f"\nRanking creado: {archivo_salida}")
    print("Genes:", len(ranking))
    print(
        "Rango mean MODZ:",
        ranking["mean_modz"].min(),
        "a",
        ranking["mean_modz"].max()
    )

# ---------------------------------------------------------------------
# 5. Guardar resumen
# ---------------------------------------------------------------------

resumen = pd.DataFrame(resumen)

resumen.to_csv(
    "RANKINGS_resumen.csv",
    index=False
)

print("\n" + "=" * 80)
print("ARCHIVOS GENERADOS")
print("=" * 80)

print("- RANKING_estradiol.csv")
print("- RANKING_calcitriol.csv")
print("- RANKING_calcipotriol.csv")
print("- RANKINGS_resumen.csv")

print("\nFIN")