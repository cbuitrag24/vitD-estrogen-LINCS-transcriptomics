import pandas as pd

archivo = "ENRIQUECIMIENTO_solapamiento_resultados.csv"

x = pd.read_csv(archivo)

significativos = x[x["q_BH"] < 0.05].copy()

columnas = [
    "analysis_space",
    "comparison",
    "direction",
    "observed_overlap",
    "fold_enrichment",
    "q_BH",
]

print("\nRESULTADOS SIGNIFICATIVOS DEL PASO 23")
print("=" * 80)

print(
    significativos[columnas]
    .sort_values(["analysis_space", "comparison", "direction"])
    .to_string(index=False)
)

print("\nNúmero total de resultados significativos:", len(significativos))

significativos[columnas].to_csv(
    "RESULTADOS_SIGNIFICATIVOS_PASO23.csv",
    index=False
)