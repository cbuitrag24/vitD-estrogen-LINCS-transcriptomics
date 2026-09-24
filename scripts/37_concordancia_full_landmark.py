import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr

print("=" * 88)
print("CONCORDANCIA GLOBAL GSEA FULL vs LANDMARK")
print("=" * 88)

# ================================================================
# 1. CARGAR COMPARACIÓN YA GENERADA
# ================================================================

df = pd.read_csv(
    "GSEA_FULL_vs_LANDMARK_comparacion.csv"
)

# Sólo Hallmarks que pudieron analizarse en Landmark
df = df[
    df["analizado_landmark"] == True
].copy()

resultados = []

# ================================================================
# 2. ANÁLISIS POR TRATAMIENTO
# ================================================================

for tratamiento in sorted(df["tratamiento"].unique()):

    sub = df[
        df["tratamiento"] == tratamiento
    ].copy()

    n = len(sub)

    # Correlaciones de NES
    rho, p_spearman = spearmanr(
        sub["NES_full"],
        sub["NES_landmark"]
    )

    r, p_pearson = pearsonr(
        sub["NES_full"],
        sub["NES_landmark"]
    )

    # Concordancia direccional
    misma_direccion = (
        np.sign(sub["NES_full"])
        ==
        np.sign(sub["NES_landmark"])
    )

    n_misma = int(misma_direccion.sum())
    porcentaje_misma = 100 * n_misma / n

    # Hallmarks significativos en FULL que sí eran evaluables
    full_sig = sub[
        sub["q_full"] < 0.05
    ]

    n_full_sig = len(full_sig)

    # De esos, cuántos siguen significativos y misma dirección
    retenidos = full_sig[
        (full_sig["q_landmark"] < 0.05)
        &
        (
            np.sign(full_sig["NES_full"])
            ==
            np.sign(full_sig["NES_landmark"])
        )
    ]

    n_retenidos = len(retenidos)

    if n_full_sig > 0:
        porcentaje_retencion = (
            100 * n_retenidos / n_full_sig
        )
    else:
        porcentaje_retencion = np.nan

    # Significativos sólo en Landmark
    solo_landmark = sub[
        (sub["q_full"] >= 0.05)
        &
        (sub["q_landmark"] < 0.05)
    ]

    resultados.append({
        "tratamiento": tratamiento,
        "hallmarks_evaluables": n,
        "spearman_rho": rho,
        "spearman_p": p_spearman,
        "pearson_r": r,
        "pearson_p": p_pearson,
        "misma_direccion_n": n_misma,
        "misma_direccion_porcentaje": porcentaje_misma,
        "full_significativos_evaluables": n_full_sig,
        "retenidos_landmark": n_retenidos,
        "retencion_porcentaje": porcentaje_retencion,
        "solo_landmark": len(solo_landmark)
    })

# ================================================================
# 3. GUARDAR Y MOSTRAR
# ================================================================

resumen = pd.DataFrame(resultados)

resumen.to_csv(
    "GSEA_FULL_vs_LANDMARK_concordancia_global.csv",
    index=False
)

print("\nRESULTADOS:")
print(
    resumen.to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}"
    )
)

print("\nArchivo generado:")
print(
    "- GSEA_FULL_vs_LANDMARK_concordancia_global.csv"
)

print("\nFIN")