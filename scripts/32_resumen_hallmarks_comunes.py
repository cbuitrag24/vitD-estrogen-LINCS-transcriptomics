import pandas as pd

print("=" * 80)
print("RESUMEN EXACTO DE HALLMARKS COMUNES A LOS TRES TRATAMIENTOS")
print("=" * 80)

archivo = "HALLMARK_comunes_3_tratamientos.csv"

df = pd.read_csv(archivo)

columnas = [
    "hallmark",
    "NES_estradiol",
    "q_estradiol",
    "NES_calcitriol",
    "q_calcitriol",
    "NES_calcipotriol",
    "q_calcipotriol",
    "patron"
]

faltantes = [c for c in columnas if c not in df.columns]

if faltantes:
    raise ValueError(
        "Faltan columnas esperadas: " + ", ".join(faltantes)
    )

tabla = df[columnas].copy()

tabla = tabla.sort_values(
    "hallmark"
).reset_index(drop=True)

tabla.to_csv(
    "HALLMARK_comunes_3_resumen_exacto.csv",
    index=False
)

print("\nNúmero de Hallmarks comunes:", len(tabla))

print("\nResultados exactos:\n")

print(
    tabla.to_string(
        index=False,
        float_format=lambda x: f"{x:.6g}"
    )
)

print("\nArchivo generado:")
print("- HALLMARK_comunes_3_resumen_exacto.csv")

print("\nFIN")