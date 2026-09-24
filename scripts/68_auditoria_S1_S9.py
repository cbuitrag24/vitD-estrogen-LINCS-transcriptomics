import pandas as pd
from pathlib import Path

# ============================================================
# AUDITORÍA FINAL — Supplementary Tables S1–S9
# No modifica ningún archivo.
# ============================================================

BASE = Path(".")

expected_rows = {
    "Supplementary_Table_S1.csv": 12,

    "Supplementary_Table_S2A.csv": 18,
    "Supplementary_Table_S2B.csv": 12,

    "Supplementary_Table_S3.csv": 4408,

    "Supplementary_Table_S4.csv": 150,
    "Supplementary_Table_S5.csv": 105,

    "Supplementary_Table_S6A.csv": 18,
    "Supplementary_Table_S6B.csv": 168,
    "Supplementary_Table_S6C.csv": 46,

    "Supplementary_Table_S7A.csv": 39,
    "Supplementary_Table_S7B.csv": 13,

    "Supplementary_Table_S8A.csv": 9,
    "Supplementary_Table_S8B.csv": 27,
    "Supplementary_Table_S8C.csv": 36,

    "Supplementary_Table_S9A.csv": 7,
    "Supplementary_Table_S9B.csv": 7,
    "Supplementary_Table_S9C.csv": 10,
}

print("=" * 72)
print("AUDITORÍA FINAL DE SUPPLEMENTARY TABLES S1–S9")
print("=" * 72)

all_ok = True

# ============================================================
# 1. Existencia y número de filas
# ============================================================

print("\n[1] EXISTENCIA Y FILAS\n")

for filename, expected in expected_rows.items():

    path = BASE / filename

    if not path.exists():
        print(f"ERROR   {filename} — NO EXISTE")
        all_ok = False
        continue

    try:
        df = pd.read_csv(path)
        observed = len(df)

        if observed == expected:
            status = "OK"
        else:
            status = "ERROR"
            all_ok = False

        print(
            f"{status:7} {filename:40} "
            f"filas={observed:5}  esperado={expected:5}"
        )

    except Exception as e:
        print(f"ERROR   {filename} — {e}")
        all_ok = False


# ============================================================
# 2. Columnas
# ============================================================

print("\n[2] COLUMNAS\n")

for filename in expected_rows:

    path = BASE / filename

    if not path.exists():
        continue

    try:
        df = pd.read_csv(path)

        print(f"\n{filename}")
        print("  " + ", ".join(df.columns))

    except Exception as e:
        print(f"ERROR leyendo {filename}: {e}")
        all_ok = False


# ============================================================
# 3. Validaciones específicas conocidas
# ============================================================

print("\n" + "=" * 72)
print("[3] VALIDACIONES ESPECÍFICAS")
print("=" * 72)


# ------------------------------------------------------------
# S1
# ------------------------------------------------------------

df = pd.read_csv("Supplementary_Table_S1.csv")

if len(df) == 12:
    print("OK      S1: 12 filas")
else:
    print("ERROR   S1: número de filas inesperado")
    all_ok = False


# ------------------------------------------------------------
# S2
# ------------------------------------------------------------

s2a = pd.read_csv("Supplementary_Table_S2A.csv")
s2b = pd.read_csv("Supplementary_Table_S2B.csv")

if len(s2a) == 18 and len(s2b) == 12:
    print("OK      S2: 18 filas descriptivas + 12 comparaciones exactas")
else:
    print("ERROR   S2: estructura inesperada")
    all_ok = False


# ------------------------------------------------------------
# S3
# ------------------------------------------------------------

s3 = pd.read_csv("Supplementary_Table_S3.csv")

if len(s3) == 4408:
    print("OK      S3: 4408 filas")

    if "Treatment" in s3.columns and "Direction" in s3.columns:
        counts = (
            s3.groupby(["Treatment", "Direction"])
            .size()
            .to_dict()
        )

        expected_s3 = {
            ("Calcipotriol", "Down"): 897,
            ("Calcipotriol", "Up"): 1063,
            ("Calcitriol", "Down"): 716,
            ("Calcitriol", "Up"): 548,
            ("Estradiol", "Down"): 402,
            ("Estradiol", "Up"): 782,
        }

        if counts == expected_s3:
            print("OK      S3: conteos 4/4 verificados")
        else:
            print("ERROR   S3: conteos 4/4 no coinciden")
            print(counts)
            all_ok = False

else:
    print("ERROR   S3: número de filas inesperado")
    all_ok = False


# ------------------------------------------------------------
# S4
# ------------------------------------------------------------

s4 = pd.read_csv("Supplementary_Table_S4.csv")

if len(s4) == 150:
    print("OK      S4: 150 filas = 50 Hallmarks × 3 tratamientos")
else:
    print("ERROR   S4: estructura inesperada")
    all_ok = False


# ------------------------------------------------------------
# S5
# ------------------------------------------------------------

s5 = pd.read_csv("Supplementary_Table_S5.csv")

if len(s5) == 105:
    print("OK      S5: 105 filas = 35 Hallmarks × 3 tratamientos")
else:
    print("ERROR   S5: estructura inesperada")
    all_ok = False


# ------------------------------------------------------------
# S6
# ------------------------------------------------------------

s6a = pd.read_csv("Supplementary_Table_S6A.csv")
s6b = pd.read_csv("Supplementary_Table_S6B.csv")
s6c = pd.read_csv("Supplementary_Table_S6C.csv")

if (
    len(s6a) == 18
    and len(s6b) == 168
    and len(s6c) == 46
):
    print("OK      S6: 18 + 168 + 46 filas")
else:
    print("ERROR   S6: estructura inesperada")
    all_ok = False


# ------------------------------------------------------------
# S7
# ------------------------------------------------------------

s7a = pd.read_csv("Supplementary_Table_S7A.csv")
s7b = pd.read_csv("Supplementary_Table_S7B.csv")

if len(s7a) == 39 and len(s7b) == 13:
    print("OK      S7: 39 + 13 filas")

    if "Context_robust_all_3_treatments" in s7b.columns:
        robust_n = int(
            s7b["Context_robust_all_3_treatments"].sum()
        )

        if robust_n == 8:
            print("OK      S7: 8 genes contextualmente robustos")
        else:
            print(
                f"ERROR   S7: se encontraron {robust_n} "
                "genes robustos; esperado=8"
            )
            all_ok = False
else:
    print("ERROR   S7: estructura inesperada")
    all_ok = False


# ------------------------------------------------------------
# S8
# ------------------------------------------------------------

s8a = pd.read_csv("Supplementary_Table_S8A.csv")
s8b = pd.read_csv("Supplementary_Table_S8B.csv")
s8c = pd.read_csv("Supplementary_Table_S8C.csv")

if (
    len(s8a) == 9
    and len(s8b) == 27
    and len(s8c) == 36
):
    print("OK      S8: 9 + 27 + 36 filas")

    if "L1000_measurement_type" in s8a.columns:

        inferred = (
            s8a["L1000_measurement_type"]
            == "Inferred"
        ).sum()

        direct = (
            s8a["L1000_measurement_type"]
            == "Landmark_direct"
        ).sum()

        if inferred == 8 and direct == 1:
            print(
                "OK      S8: 8 genes inferidos + "
                "1 Landmark directo"
            )
        else:
            print(
                "ERROR   S8: tipo de medición inesperado"
            )
            all_ok = False

else:
    print("ERROR   S8: estructura inesperada")
    all_ok = False


# ------------------------------------------------------------
# S9
# ------------------------------------------------------------

s9a = pd.read_csv("Supplementary_Table_S9A.csv")
s9b = pd.read_csv("Supplementary_Table_S9B.csv")
s9c = pd.read_csv("Supplementary_Table_S9C.csv")

if (
    len(s9a) == 7
    and len(s9b) == 7
    and len(s9c) == 10
):
    print("OK      S9: 7 + 7 + 10 filas")

    if "Cell_line" in s9a.columns:

        cell_counts = (
            s9a.groupby("Cell_line")
            .size()
            .to_dict()
        )

        expected_cells = {
            "A549": 2,
            "PC3": 4,
            "VCAP": 1,
        }

        if cell_counts == expected_cells:
            print(
                "OK      S9: 7 bloques pareados "
                "(A549=2, PC3=4, VCAP=1)"
            )
        else:
            print(
                "ERROR   S9: distribución de bloques inesperada"
            )
            print(cell_counts)
            all_ok = False

else:
    print("ERROR   S9: estructura inesperada")
    all_ok = False


# ============================================================
# 4. Resumen final
# ============================================================

print("\n" + "=" * 72)
print("RESULTADO FINAL DE LA AUDITORÍA")
print("=" * 72)

if all_ok:
    print("\n✓ TODAS LAS VALIDACIONES PASARON CORRECTAMENTE.")
    print("✓ Los 17 archivos S1–S9 están presentes.")
    print("✓ Los conteos de filas coinciden con los resultados verificados.")
    print("✓ Las validaciones específicas pasaron.")
    print("\nNo se modificó ningún archivo.")
else:
    print("\n✗ SE DETECTARON UNA O MÁS INCONSISTENCIAS.")
    print("Revisar los mensajes marcados como ERROR.")

print()