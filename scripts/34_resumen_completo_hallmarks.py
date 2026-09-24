import pandas as pd

print("=" * 88)
print("RESUMEN COMPLETO DE LOS 50 HALLMARKS")
print("=" * 88)

archivo = "HALLMARK_comparacion_tratamientos.csv"

df = pd.read_csv(archivo)

columnas_necesarias = [
    "hallmark",
    "NES_estradiol",
    "q_estradiol",
    "sig_estradiol",
    "NES_calcitriol",
    "q_calcitriol",
    "sig_calcitriol",
    "NES_calcipotriol",
    "q_calcipotriol",
    "sig_calcipotriol",
    "n_tratamientos_significativos",
    "patron"
]

faltantes = [c for c in columnas_necesarias if c not in df.columns]

if faltantes:
    raise ValueError(
        "Faltan columnas esperadas: " + ", ".join(faltantes)
    )

df = df[columnas_necesarias].copy()

# Ordenar primero por patrón y luego alfabéticamente
orden_patrones = [
    "comun_3_misma_direccion",
    "comun_3_direccion_divergente",
    "estradiol_calcitriol_misma_direccion",
    "estradiol_calcitriol_opuesta",
    "estradiol_calcipotriol_misma_direccion",
    "estradiol_calcipotriol_opuesta",
    "calcitriol_calcipotriol_misma_direccion",
    "calcitriol_calcipotriol_opuesta",
    "solo_estradiol",
    "solo_calcitriol",
    "solo_calcipotriol",
    "no_significativo"
]

df["patron"] = pd.Categorical(
    df["patron"],
    categories=orden_patrones,
    ordered=True
)

df = df.sort_values(
    ["patron", "hallmark"]
).reset_index(drop=True)

# Guardar tabla ordenada
df.to_csv(
    "HALLMARK_resumen_completo_50.csv",
    index=False
)

print("\nNúmero total de Hallmarks:", len(df))

# Mostrar cada grupo por separado
for patron in orden_patrones:

    subset = df[df["patron"] == patron].copy()

    if len(subset) == 0:
        continue

    print("\n" + "=" * 88)
    print(f"PATRÓN: {patron}")
    print(f"N = {len(subset)}")
    print("=" * 88)

    for _, fila in subset.iterrows():

        print(f"\n{fila['hallmark']}")

        print(
            f"  Estradiol:    NES = {fila['NES_estradiol']:.6f}   "
            f"q = {fila['q_estradiol']:.8g}   "
            f"sig = {fila['sig_estradiol']}"
        )

        print(
            f"  Calcitriol:   NES = {fila['NES_calcitriol']:.6f}   "
            f"q = {fila['q_calcitriol']:.8g}   "
            f"sig = {fila['sig_calcitriol']}"
        )

        print(
            f"  Calcipotriol: NES = {fila['NES_calcipotriol']:.6f}   "
            f"q = {fila['q_calcipotriol']:.8g}   "
            f"sig = {fila['sig_calcipotriol']}"
        )

print("\n" + "=" * 88)
print("ARCHIVO GENERADO")
print("=" * 88)
print("- HALLMARK_resumen_completo_50.csv")

print("\nFIN")