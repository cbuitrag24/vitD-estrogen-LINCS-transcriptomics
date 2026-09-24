import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Supplementary Table S1 — publication-ready editorial version
# The audited source CSV is read only and is never modified.
# ------------------------------------------------------------

INPUT_FILE = Path("Supplementary_Table_S1.csv")
OUTPUT_FILE = Path("Supplementary_Table_S1_PUBLICATION.csv")

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"No se encontró el archivo fuente: {INPUT_FILE.resolve()}"
    )

# Read audited source data
df = pd.read_csv(INPUT_FILE)

# Verify expected structure before doing anything
expected_columns = [
    "Cell_line",
    "Experiment_prefix",
    "Matched_block",
    "Perturbagen",
    "LINCS_pert_id",
    "LINCS_sig_id",
    "Dose",
    "Dose_unit",
    "Exposure_time",
    "Time_unit",
    "N_constituent_profiles",
    "Constituent_profile_IDs",
]

missing = [col for col in expected_columns if col not in df.columns]

if missing:
    raise ValueError(
        "Faltan columnas esperadas en S1: " + ", ".join(missing)
    )

if len(df) != 12:
    raise ValueError(
        f"S1 debería contener 12 filas, pero contiene {len(df)}."
    )

# Work on a copy — never modify the audited dataframe
pub = df.copy()

# Publication-friendly column headings
pub = pub.rename(
    columns={
        "Cell_line": "Cell line",
        "Experiment_prefix": "Experimental prefix",
        "Matched_block": "Matched experimental block",
        "Perturbagen": "Perturbagen",
        "LINCS_pert_id": "LINCS perturbagen ID",
        "LINCS_sig_id": "LINCS signature ID",
        "Dose": "Dose",
        "Dose_unit": "Dose unit",
        "Exposure_time": "Exposure time",
        "Time_unit": "Time unit",
        "N_constituent_profiles": "N constituent profiles",
        "Constituent_profile_IDs": "Constituent profile IDs",
    }
)

# Standardise display terminology only.
# Numerical values and LINCS identifiers remain unchanged.
pub["Perturbagen"] = pub["Perturbagen"].replace(
    {
        "estradiol": "Estradiol",
        "calcitriol": "Calcitriol",
        "calcipotriol": "Calcipotriol",
    }
)

# Keep the audited row order and all source identifiers
pub.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

print()
print("S1 PUBLICATION VERSION CREATED CORRECTLY")
print("----------------------------------------")
print(f"Source file: {INPUT_FILE}")
print(f"Output file: {OUTPUT_FILE}")
print(f"Rows: {len(pub)}")
print(f"Columns: {len(pub.columns)}")
print()
print("The audited source CSV was NOT modified.")