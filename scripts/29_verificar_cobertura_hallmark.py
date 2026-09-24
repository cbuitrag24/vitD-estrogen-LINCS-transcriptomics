import pandas as pd

print("=" * 80)
print("VERIFICACIÓN DE COBERTURA HALLMARK")
print("=" * 80)

# Genes presentes en nuestro dataset
genes = pd.read_csv("dataset_primario_genes.csv")

genes_lincs = set(
    genes["pr_gene_symbol"]
    .dropna()
    .astype(str)
    .str.strip()
)

print("\nGenes LINCS únicos:", len(genes_lincs))

# Leer Hallmark GMT
hallmarks = {}

with open(
    "h.all.v2026.1.Hs.symbols.gmt",
    "r",
    encoding="utf-8"
) as f:

    for linea in f:

        partes = linea.strip().split("\t")

        if len(partes) < 3:
            continue

        nombre = partes[0]
        genes_set = set(partes[2:])

        hallmarks[nombre] = genes_set

print("Hallmark gene sets:", len(hallmarks))

# Comprobar cobertura
resultados = []

for nombre, genes_set in hallmarks.items():

    presentes = genes_set & genes_lincs

    resultados.append({
        "hallmark": nombre,
        "genes_hallmark": len(genes_set),
        "genes_en_LINCS": len(presentes),
        "cobertura_porcentaje":
            100 * len(presentes) / len(genes_set)
    })

resultado = pd.DataFrame(resultados)

resultado.to_csv(
    "HALLMARK_cobertura_LINCS.csv",
    index=False
)

print("\nResumen de cobertura:")
print(
    resultado[
        ["genes_en_LINCS", "cobertura_porcentaje"]
    ].describe()
)

print("\nGene sets con menor cobertura:")

print(
    resultado
    .sort_values("cobertura_porcentaje")
    .head(10)
    .to_string(index=False)
)

print("\nNúmero de Hallmarks con al menos 15 genes presentes:")
print((resultado["genes_en_LINCS"] >= 15).sum())

print("\nArchivo generado:")
print("- HALLMARK_cobertura_LINCS.csv")

print("\nFIN")