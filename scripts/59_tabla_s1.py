import pandas as pd

# ============================================================
# Supplementary Table S1
# Detailed metadata of the 12 primary LINCS L1000 signatures
# ============================================================

INPUT_FILE = "dataset_primario_metadata.csv"
OUTPUT_FILE = "Supplementary_Table_S1.csv"

# Leer los metadatos primarios ya validados
df = pd.read_csv(INPUT_FILE)

# Seleccionar únicamente información necesaria para la
# trazabilidad del dataset primario
s1 = df[
    [
        "cell_id",
        "experiment_prefix",
        "block_id",
        "pert_iname",
        "pert_id",
        "sig_id",
        "pert_dose",
        "pert_dose_unit",
        "pert_time",
        "pert_time_unit",
        "n_replicas",
        "distil_id",
    ]
].copy()

# Nombres destinados a la tabla suplementaria
s1 = s1.rename(
    columns={
        "cell_id": "Cell_line",
        "experiment_prefix": "Experiment_prefix",
        "block_id": "Matched_block",
        "pert_iname": "Perturbagen",
        "pert_id": "LINCS_pert_id",
        "sig_id": "LINCS_sig_id",
        "pert_dose": "Dose",
        "pert_dose_unit": "Dose_unit",
        "pert_time": "Exposure_time",
        "pert_time_unit": "Time_unit",
        "n_replicas": "N_constituent_profiles",
        "distil_id": "Constituent_profile_IDs",
    }
)

# Ordenar para facilitar la lectura de los cuatro bloques
treatment_order = {
    "estradiol": 0,
    "calcitriol": 1,
    "calcipotriol": 2,
}

s1["_treatment_order"] = s1["Perturbagen"].map(treatment_order)

s1 = (
    s1.sort_values(
        ["Cell_line", "_treatment_order"],
        kind="stable"
    )
    .drop(columns="_treatment_order")
    .reset_index(drop=True)
)

# Controles antes de guardar
assert len(s1) == 12
assert s1["Matched_block"].nunique() == 4
assert set(s1["Perturbagen"]) == {
    "estradiol",
    "calcitriol",
    "calcipotriol",
}
assert s1["N_constituent_profiles"].sum() == 57

# Guardar
s1.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

print("\nSupplementary Table S1 generated successfully.")
print(f"Rows: {len(s1)}")
print(f"Matched blocks: {s1['Matched_block'].nunique()}")
print(f"Total constituent profiles: {s1['N_constituent_profiles'].sum()}")
print(f"Output file: {OUTPUT_FILE}\n")

print(s1.to_string(index=False))