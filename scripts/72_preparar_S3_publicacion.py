import pandas as pd

INPUT_FILE = "Supplementary_Table_S3.csv"
OUTPUT_FILE = "Supplementary_Data_Table_S3_PUBLICATION.csv"

# Read the audited source table without modifying it
df = pd.read_csv(INPUT_FILE)

# Safety check: S3 must contain exactly 4,408 rows
if len(df) != 4408:
    raise ValueError(
        f"Unexpected number of rows in S3: {len(df)}. Expected 4408."
    )

# Publication-friendly column headings
rename_columns = {
    "Gene_ID": "Gene ID",
    "Gene_symbol": "Gene symbol",
    "Gene_title": "Gene title",
    "Landmark_status": "Landmark status",
    "Treatment": "Treatment",
    "N_contexts": "N contexts",
    "N_positive": "N positive contexts",
    "N_negative": "N negative contexts",
    "Direction": "Direction",
    "Consistency": "Consistency",
    "Mean_MODZ": "Mean MODZ",
    "Median_MODZ": "Median MODZ",
    "Mean_abs_MODZ": "Mean absolute MODZ",
}

publication = df.rename(
    columns={
        key: value
        for key, value in rename_columns.items()
        if key in df.columns
    }
).copy()

# Convert Landmark status to a reader-friendly label in the
# publication copy only. The audited source remains unchanged.
if "Landmark status" in publication.columns:
    publication["Landmark status"] = publication["Landmark status"].map(
        {1: "Yes", 0: "No", "1": "Yes", "0": "No"}
    ).fillna(publication["Landmark status"])

# Save a separate publication copy
publication.to_csv(OUTPUT_FILE, index=False)

print("S3 PUBLICATION VERSION CREATED CORRECTLY")
print("----------------------------------------")
print(f"Output file: {OUTPUT_FILE}")
print(f"Rows: {len(publication)}")
print(f"Columns: {len(publication.columns)}")
print()
print("The audited Supplementary_Table_S3.csv file was NOT modified.")