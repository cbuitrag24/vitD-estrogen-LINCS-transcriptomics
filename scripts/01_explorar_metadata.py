import pandas as pd

sig_file = "data/GSE92742_Broad_LINCS_sig_info.txt.gz"
pert_file = "data/GSE92742_Broad_LINCS_pert_info.txt.gz"
gene_file = "data/GSE92742_Broad_LINCS_gene_info.txt.gz"

sig = pd.read_csv(sig_file, sep="\t", compression="gzip", low_memory=False)
pert = pd.read_csv(pert_file, sep="\t", compression="gzip", low_memory=False)
genes = pd.read_csv(gene_file, sep="\t", compression="gzip", low_memory=False)

print("SIG_INFO")
print("Dimensiones:", sig.shape)
print("Columnas:", sig.columns.tolist())

print("\nPERT_INFO")
print("Dimensiones:", pert.shape)
print("Columnas:", pert.columns.tolist())

print("\nGENE_INFO")
print("Dimensiones:", genes.shape)
print("Columnas:", genes.columns.tolist())