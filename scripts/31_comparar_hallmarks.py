import pandas as pd
import numpy as np

print("=" * 80)
print("COMPARACIÓN DE HALLMARKS ENTRE TRATAMIENTOS")
print("=" * 80)

archivo = "GSEA_HALLMARK_todos_tratamientos.csv"

df = pd.read_csv(archivo)

columnas_necesarias = {
    "tratamiento",
    "hallmark",
    "NES",
    "p_empirico",
    "q_BH"
}

if not columnas_necesarias.issubset(df.columns):
    raise ValueError(
        "ERROR: faltan columnas necesarias en "
        "GSEA_HALLMARK_todos_tratamientos.csv"
    )

tratamientos = [
    "estradiol",
    "calcitriol",
    "calcipotriol"
]

# ------------------------------------------------------------
# 1. Crear tabla ancha con NES y q_BH
# ------------------------------------------------------------

nes = df.pivot(
    index="hallmark",
    columns="tratamiento",
    values="NES"
)

q = df.pivot(
    index="hallmark",
    columns="tratamiento",
    values="q_BH"
)

comparacion = pd.DataFrame(index=nes.index)

for tratamiento in tratamientos:
    comparacion[f"NES_{tratamiento}"] = nes[tratamiento]
    comparacion[f"q_{tratamiento}"] = q[tratamiento]

    comparacion[f"sig_{tratamiento}"] = (
        q[tratamiento] < 0.05
    )

comparacion = comparacion.reset_index()

# ------------------------------------------------------------
# 2. Número de tratamientos significativos
# ------------------------------------------------------------

comparacion["n_tratamientos_significativos"] = (
    comparacion[
        [
            "sig_estradiol",
            "sig_calcitriol",
            "sig_calcipotriol"
        ]
    ]
    .sum(axis=1)
)

# ------------------------------------------------------------
# 3. Dirección del NES
# ------------------------------------------------------------

def direccion(valor):
    if valor > 0:
        return "positivo"
    elif valor < 0:
        return "negativo"
    else:
        return "cero"


for tratamiento in tratamientos:
    comparacion[f"direccion_{tratamiento}"] = (
        comparacion[f"NES_{tratamiento}"]
        .apply(direccion)
    )

# ------------------------------------------------------------
# 4. Clasificación de patrón
# ------------------------------------------------------------

def clasificar_fila(fila):

    sig_e = fila["sig_estradiol"]
    sig_c = fila["sig_calcitriol"]
    sig_cp = fila["sig_calcipotriol"]

    nes_e = fila["NES_estradiol"]
    nes_c = fila["NES_calcitriol"]
    nes_cp = fila["NES_calcipotriol"]

    # Tres tratamientos significativos
    if sig_e and sig_c and sig_cp:

        signos = np.sign([nes_e, nes_c, nes_cp])

        if len(set(signos)) == 1:
            return "comun_3_misma_direccion"
        else:
            return "comun_3_direccion_divergente"

    # Pares
    if sig_e and sig_c and not sig_cp:

        if np.sign(nes_e) == np.sign(nes_c):
            return "estradiol_calcitriol_misma_direccion"
        else:
            return "estradiol_calcitriol_opuesta"

    if sig_e and sig_cp and not sig_c:

        if np.sign(nes_e) == np.sign(nes_cp):
            return "estradiol_calcipotriol_misma_direccion"
        else:
            return "estradiol_calcipotriol_opuesta"

    if sig_c and sig_cp and not sig_e:

        if np.sign(nes_c) == np.sign(nes_cp):
            return "calcitriol_calcipotriol_misma_direccion"
        else:
            return "calcitriol_calcipotriol_opuesta"

    # Exclusivos
    if sig_e and not sig_c and not sig_cp:
        return "solo_estradiol"

    if sig_c and not sig_e and not sig_cp:
        return "solo_calcitriol"

    if sig_cp and not sig_e and not sig_c:
        return "solo_calcipotriol"

    return "no_significativo"


comparacion["patron"] = comparacion.apply(
    clasificar_fila,
    axis=1
)

# ------------------------------------------------------------
# 5. Guardar tabla completa
# ------------------------------------------------------------

comparacion = comparacion.sort_values(
    [
        "n_tratamientos_significativos",
        "hallmark"
    ],
    ascending=[False, True]
)

comparacion.to_csv(
    "HALLMARK_comparacion_tratamientos.csv",
    index=False
)

# ------------------------------------------------------------
# 6. Resumen de patrones
# ------------------------------------------------------------

resumen = (
    comparacion["patron"]
    .value_counts()
    .rename_axis("patron")
    .reset_index(name="n_hallmarks")
)

resumen.to_csv(
    "HALLMARK_resumen_patrones.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Hallmarks significativos en los tres
# ------------------------------------------------------------

comunes_3 = comparacion[
    comparacion["n_tratamientos_significativos"] == 3
].copy()

comunes_3.to_csv(
    "HALLMARK_comunes_3_tratamientos.csv",
    index=False
)

# ------------------------------------------------------------
# 8. Hallmarks exclusivos
# ------------------------------------------------------------

for tratamiento in tratamientos:

    exclusivos = comparacion[
        comparacion["patron"] == f"solo_{tratamiento}"
    ].copy()

    exclusivos.to_csv(
        f"HALLMARK_exclusivos_{tratamiento}.csv",
        index=False
    )

# ------------------------------------------------------------
# 9. Mostrar resumen
# ------------------------------------------------------------

print("\nResumen de patrones:")
print(resumen.to_string(index=False))

print("\nHallmarks significativos en los tres tratamientos:")
print(len(comunes_3))

if len(comunes_3) > 0:
    columnas_mostrar = [
        "hallmark",
        "NES_estradiol",
        "NES_calcitriol",
        "NES_calcipotriol",
        "patron"
    ]

    print(
        comunes_3[columnas_mostrar]
        .to_string(index=False)
    )

print("\nArchivos generados:")
print("- HALLMARK_comparacion_tratamientos.csv")
print("- HALLMARK_resumen_patrones.csv")
print("- HALLMARK_comunes_3_tratamientos.csv")
print("- HALLMARK_exclusivos_estradiol.csv")
print("- HALLMARK_exclusivos_calcitriol.csv")
print("- HALLMARK_exclusivos_calcipotriol.csv")

print("\nFIN")