import pandas as pd
import numpy as np

print("=" * 88)
print("COMPARACIÓN GSEA FULL vs LANDMARK")
print("=" * 88)

# ================================================================
# 1. CARGAR RESULTADOS
# ================================================================

full = pd.read_csv(
    "GSEA_HALLMARK_todos_tratamientos.csv"
)

landmark = pd.read_csv(
    "GSEA_LANDMARK_todos_tratamientos.csv"
)

columnas_full = {
    "tratamiento",
    "hallmark",
    "NES",
    "q_BH"
}

columnas_landmark = {
    "tratamiento",
    "hallmark",
    "NES",
    "q_BH"
}

if not columnas_full.issubset(full.columns):
    raise ValueError(
        "Faltan columnas esperadas en resultados FULL"
    )

if not columnas_landmark.issubset(landmark.columns):
    raise ValueError(
        "Faltan columnas esperadas en resultados LANDMARK"
    )

# ================================================================
# 2. RENOMBRAR COLUMNAS
# ================================================================

full = full[
    ["tratamiento", "hallmark", "NES", "q_BH"]
].copy()

landmark = landmark[
    ["tratamiento", "hallmark", "NES", "q_BH"]
].copy()

full = full.rename(
    columns={
        "NES": "NES_full",
        "q_BH": "q_full"
    }
)

landmark = landmark.rename(
    columns={
        "NES": "NES_landmark",
        "q_BH": "q_landmark"
    }
)

# ================================================================
# 3. UNIR
# ================================================================

comparacion = full.merge(
    landmark,
    on=["tratamiento", "hallmark"],
    how="left"
)

comparacion["analizado_landmark"] = (
    comparacion["NES_landmark"].notna()
)

comparacion["sig_full"] = (
    comparacion["q_full"] < 0.05
)

comparacion["sig_landmark"] = (
    comparacion["q_landmark"] < 0.05
)

# ================================================================
# 4. DIRECCIONES
# ================================================================

def signo(valor):

    if pd.isna(valor):
        return np.nan

    if valor > 0:
        return 1

    if valor < 0:
        return -1

    return 0


comparacion["signo_full"] = (
    comparacion["NES_full"].apply(signo)
)

comparacion["signo_landmark"] = (
    comparacion["NES_landmark"].apply(signo)
)

comparacion["misma_direccion"] = (
    comparacion["signo_full"]
    ==
    comparacion["signo_landmark"]
).astype("boolean")

comparacion.loc[
    ~comparacion["analizado_landmark"],
    "misma_direccion"
] = pd.NA

# ================================================================
# 5. CLASIFICACIÓN DE ROBUSTEZ
# ================================================================

def clasificar(fila):

    if not fila["analizado_landmark"]:
        return "no_evaluable_landmark"

    if (
        fila["sig_full"]
        and fila["sig_landmark"]
        and fila["misma_direccion"]
    ):
        return "robusto_full_landmark"

    if (
        fila["sig_full"]
        and fila["sig_landmark"]
        and not fila["misma_direccion"]
    ):
        return "significativo_ambos_direccion_opuesta"

    if (
        fila["sig_full"]
        and not fila["sig_landmark"]
    ):
        return "solo_full"

    if (
        not fila["sig_full"]
        and fila["sig_landmark"]
    ):
        return "solo_landmark"

    return "no_significativo_ambos"


comparacion["robustez"] = (
    comparacion.apply(
        clasificar,
        axis=1
    )
)

# ================================================================
# 6. GUARDAR TABLA COMPLETA
# ================================================================

comparacion.to_csv(
    "GSEA_FULL_vs_LANDMARK_comparacion.csv",
    index=False
)

# ================================================================
# 7. RESUMEN POR TRATAMIENTO
# ================================================================

resumen = (
    comparacion
    .groupby(
        ["tratamiento", "robustez"]
    )
    .size()
    .reset_index(name="n_hallmarks")
)

resumen.to_csv(
    "GSEA_FULL_vs_LANDMARK_resumen.csv",
    index=False
)

print("\nResumen por tratamiento y categoría:")
print(
    resumen.to_string(index=False)
)

# ================================================================
# 8. LAS 8 HALLMARKS COMUNES DEL ANÁLISIS FULL
# ================================================================

hallmarks_comunes = [
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
]

nucleo = comparacion[
    comparacion["hallmark"].isin(
        hallmarks_comunes
    )
].copy()

nucleo = nucleo.sort_values(
    ["hallmark", "tratamiento"]
)

nucleo.to_csv(
    "GSEA_FULL_vs_LANDMARK_nucleo_8.csv",
    index=False
)

print("\n" + "=" * 88)
print("ROBUSTEZ DE LAS 8 HALLMARKS COMUNES DEL ANÁLISIS FULL")
print("=" * 88)

for hallmark in hallmarks_comunes:

    print("\n" + hallmark)

    subset = nucleo[
        nucleo["hallmark"] == hallmark
    ]

    for _, fila in subset.iterrows():

        if fila["analizado_landmark"]:

            print(
                f"  {fila['tratamiento']}: "
                f"FULL NES={fila['NES_full']:.6f}, "
                f"q={fila['q_full']:.8g} | "
                f"LANDMARK NES={fila['NES_landmark']:.6f}, "
                f"q={fila['q_landmark']:.8g} | "
                f"{fila['robustez']}"
            )

        else:

            print(
                f"  {fila['tratamiento']}: "
                f"FULL NES={fila['NES_full']:.6f}, "
                f"q={fila['q_full']:.8g} | "
                f"LANDMARK: no evaluable (<15 genes)"
            )

# ================================================================
# 9. RESUMEN DEL NÚCLEO DE 8
# ================================================================

resumen_nucleo = (
    nucleo
    .groupby("hallmark")
    .agg(
        tratamientos_robustos=(
            "robustez",
            lambda x:
                (x == "robusto_full_landmark").sum()
        ),
        tratamientos_no_evaluables=(
            "robustez",
            lambda x:
                (x == "no_evaluable_landmark").sum()
        ),
        tratamientos_solo_full=(
            "robustez",
            lambda x:
                (x == "solo_full").sum()
        ),
        tratamientos_solo_landmark=(
            "robustez",
            lambda x:
                (x == "solo_landmark").sum()
        )
    )
    .reset_index()
)

resumen_nucleo.to_csv(
    "GSEA_FULL_vs_LANDMARK_nucleo_8_resumen.csv",
    index=False
)

print("\n" + "=" * 88)
print("RESUMEN DEL NÚCLEO DE 8")
print("=" * 88)

print(
    resumen_nucleo.to_string(index=False)
)

print("\nArchivos generados:")
print("- GSEA_FULL_vs_LANDMARK_comparacion.csv")
print("- GSEA_FULL_vs_LANDMARK_resumen.csv")
print("- GSEA_FULL_vs_LANDMARK_nucleo_8.csv")
print("- GSEA_FULL_vs_LANDMARK_nucleo_8_resumen.csv")

print("\nFIN")