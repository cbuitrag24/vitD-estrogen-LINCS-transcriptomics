import pandas as pd

archivo = "HALLMARK_comunes_3_resumen_exacto.csv"

df = pd.read_csv(archivo)

print("\n" + "=" * 72)
print("8 HALLMARKS COMUNES A ESTRADIOL, CALCITRIOL Y CALCIPOTRIOL")
print("=" * 72)

for _, fila in df.iterrows():

    print("\n" + fila["hallmark"])

    print(
        f"  Estradiol:    NES = {fila['NES_estradiol']:.6f}   "
        f"q = {fila['q_estradiol']:.8g}"
    )

    print(
        f"  Calcitriol:   NES = {fila['NES_calcitriol']:.6f}   "
        f"q = {fila['q_calcitriol']:.8g}"
    )

    print(
        f"  Calcipotriol: NES = {fila['NES_calcipotriol']:.6f}   "
        f"q = {fila['q_calcipotriol']:.8g}"
    )

print("\n" + "=" * 72)
print("FIN")