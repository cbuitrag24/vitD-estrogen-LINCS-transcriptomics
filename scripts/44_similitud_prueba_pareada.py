import itertools
import numpy as np
import pandas as pd

print("=" * 100)
print("SIMILITUD TRANSCRIPCIONAL - COMPARACIONES PAREADAS EXACTAS")
print("=" * 100)

# ================================================================
# 1. ARCHIVO DE ENTRADA
# ================================================================

ARCHIVO = "SIMILITUD_robustez_full_vs_landmark.csv"

df = pd.read_csv(ARCHIVO)

print("\nFilas:", len(df))
print("Contextos:", sorted(df["cell_id"].unique()))
print("Comparaciones:", sorted(df["comparison"].unique()))

# ================================================================
# 2. DEFINICIONES
# ================================================================

ANALISIS = {
    "full": {
        "pearson": "pearson_full",
        "spearman": "spearman_full",
        "cosine": "cosine_full"
    },
    "landmark": {
        "pearson": "pearson_landmark",
        "spearman": "spearman_landmark",
        "cosine": "cosine_landmark"
    }
}

VD_PAIR = "calcitriol_vs_calcipotriol"

E2_PAIRS = [
    "estradiol_vs_calcitriol",
    "estradiol_vs_calcipotriol"
]

# ================================================================
# 3. FUNCIÓN EXACTA SIGN-FLIP
# ================================================================

def exact_sign_flip_test(differences):
    """
    Prueba exacta bilateral de sign-flip sobre la media
    de diferencias pareadas.

    H0: la distribución de las diferencias es simétrica
    alrededor de cero.

    Con n=4 hay 2^4 = 16 configuraciones exactas.
    """

    d = np.asarray(differences, dtype=float)

    observed = np.mean(d)

    permuted = []

    for signs in itertools.product([-1, 1], repeat=len(d)):
        signs = np.asarray(signs)
        permuted.append(np.mean(d * signs))

    permuted = np.asarray(permuted)

    p_two_sided = np.mean(
        np.abs(permuted) >=
        (abs(observed) - 1e-12)
    )

    return observed, p_two_sided, len(permuted)


# ================================================================
# 4. EXTRAER VALORES POR CONTEXTO
# ================================================================

resultados = []
detalles = []

for espacio, metricas in ANALISIS.items():

    for metrica, columna in metricas.items():

        pivot = df.pivot(
            index="cell_id",
            columns="comparison",
            values=columna
        )

        columnas_necesarias = [
            VD_PAIR,
            *E2_PAIRS
        ]

        faltantes = [
            c for c in columnas_necesarias
            if c not in pivot.columns
        ]

        if faltantes:
            raise ValueError(
                f"Faltan comparaciones para {espacio}/{metrica}: "
                f"{faltantes}"
            )

        # --------------------------------------------------------
        # A. VitD vs E2-calcitriol
        # --------------------------------------------------------

        for comparacion_e2 in E2_PAIRS:

            vd = pivot[VD_PAIR]
            e2 = pivot[comparacion_e2]

            diff = vd - e2

            media_diff, p_exacta, n_perm = (
                exact_sign_flip_test(diff.values)
            )

            mediana_diff = float(
                np.median(diff.values)
            )

            n_vd_mayor = int(
                np.sum(diff.values > 0)
            )

            n_vd_menor = int(
                np.sum(diff.values < 0)
            )

            resultados.append({
                "espacio": espacio,
                "metrica": metrica,
                "contraste":
                    f"{VD_PAIR}_menos_{comparacion_e2}",
                "n_bloques": len(diff),
                "media_similitud_VitD": float(vd.mean()),
                "media_similitud_comparador": float(e2.mean()),
                "diferencia_media": float(media_diff),
                "diferencia_mediana": mediana_diff,
                "bloques_VitD_mayor": n_vd_mayor,
                "bloques_VitD_menor": n_vd_menor,
                "n_permutaciones_exactas": n_perm,
                "p_exacta_bilateral": float(p_exacta)
            })

            for cell in pivot.index:

                detalles.append({
                    "espacio": espacio,
                    "metrica": metrica,
                    "contraste":
                        f"{VD_PAIR}_menos_{comparacion_e2}",
                    "cell_id": cell,
                    "similitud_VitD": float(vd.loc[cell]),
                    "similitud_comparador": float(e2.loc[cell]),
                    "diferencia_pareada": float(diff.loc[cell])
                })

        # --------------------------------------------------------
        # B. Comparación contra el promedio de las dos parejas E2
        # --------------------------------------------------------

        e2_promedio = (
            pivot[E2_PAIRS[0]]
            +
            pivot[E2_PAIRS[1]]
        ) / 2.0

        vd = pivot[VD_PAIR]

        diff = vd - e2_promedio

        media_diff, p_exacta, n_perm = (
            exact_sign_flip_test(diff.values)
        )

        resultados.append({
            "espacio": espacio,
            "metrica": metrica,
            "contraste":
                "calcitriol_vs_calcipotriol_menos_promedio_E2",
            "n_bloques": len(diff),
            "media_similitud_VitD": float(vd.mean()),
            "media_similitud_comparador": float(
                e2_promedio.mean()
            ),
            "diferencia_media": float(media_diff),
            "diferencia_mediana": float(
                np.median(diff.values)
            ),
            "bloques_VitD_mayor": int(
                np.sum(diff.values > 0)
            ),
            "bloques_VitD_menor": int(
                np.sum(diff.values < 0)
            ),
            "n_permutaciones_exactas": n_perm,
            "p_exacta_bilateral": float(p_exacta)
        })

        for cell in pivot.index:

            detalles.append({
                "espacio": espacio,
                "metrica": metrica,
                "contraste":
                    "calcitriol_vs_calcipotriol_menos_promedio_E2",
                "cell_id": cell,
                "similitud_VitD": float(vd.loc[cell]),
                "similitud_comparador": float(
                    e2_promedio.loc[cell]
                ),
                "diferencia_pareada": float(
                    diff.loc[cell]
                )
            })


resultados_df = pd.DataFrame(resultados)
detalles_df = pd.DataFrame(detalles)

# ================================================================
# 5. AJUSTE BH
# ================================================================

def bh_adjust(pvalues):

    p = np.asarray(pvalues, dtype=float)
    n = len(p)

    order = np.argsort(p)
    ranked = p[order]

    adjusted = (
        ranked * n /
        np.arange(1, n + 1)
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.minimum(
        adjusted,
        1.0
    )

    output = np.empty(n)
    output[order] = adjusted

    return output


resultados_df["q_BH_global"] = bh_adjust(
    resultados_df["p_exacta_bilateral"].values
)

# Ajuste también dentro de cada espacio:
# Full y Landmark se consideran análisis de sensibilidad relacionados.

resultados_df["q_BH_dentro_espacio"] = np.nan

for espacio in resultados_df["espacio"].unique():

    mask = (
        resultados_df["espacio"]
        ==
        espacio
    )

    resultados_df.loc[
        mask,
        "q_BH_dentro_espacio"
    ] = bh_adjust(
        resultados_df.loc[
            mask,
            "p_exacta_bilateral"
        ].values
    )

# ================================================================
# 6. GUARDAR
# ================================================================

resultados_df.to_csv(
    "SIMILITUD_prueba_pareada_exacta.csv",
    index=False
)

detalles_df.to_csv(
    "SIMILITUD_prueba_pareada_detalle.csv",
    index=False
)

# ================================================================
# 7. MOSTRAR RESULTADOS
# ================================================================

print("\n" + "=" * 100)
print("RESULTADOS EXACTOS")
print("=" * 100)

for espacio in ["full", "landmark"]:

    print("\n" + "-" * 100)
    print(espacio.upper())
    print("-" * 100)

    sub = resultados_df[
        resultados_df["espacio"] == espacio
    ]

    print(
        sub[
            [
                "metrica",
                "contraste",
                "media_similitud_VitD",
                "media_similitud_comparador",
                "diferencia_media",
                "diferencia_mediana",
                "bloques_VitD_mayor",
                "bloques_VitD_menor",
                "p_exacta_bilateral",
                "q_BH_dentro_espacio"
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}"
        )
    )

# ================================================================
# 8. DETALLE POR CONTEXTO
# ================================================================

print("\n" + "=" * 100)
print("DIFERENCIAS PAREADAS POR CONTEXTO")
print("=" * 100)

print(
    detalles_df.to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}"
    )
)

# ================================================================
# 9. NOTA METODOLÓGICA
# ================================================================

print("\n" + "=" * 100)
print("NOTA METODOLÓGICA")
print("=" * 100)

print(
    "Las comparaciones son pareadas por contexto/bloque experimental."
)

print(
    "La prueba exacta bilateral utiliza todas las configuraciones "
    "de cambio de signo de las diferencias pareadas."
)

print(
    "Con cuatro bloques existen solamente 2^4 = 16 configuraciones, "
    "por lo que la resolución mínima posible de una p bilateral exacta "
    "es limitada."
)

print(
    "Los resultados deben interpretarse junto con las diferencias "
    "por contexto y no únicamente mediante el p-value."
)

print(
    "El contraste frente al promedio de las dos comparaciones con "
    "estradiol es un resumen descriptivo complementario."
)

print(
    "Full y Landmark no son datasets independientes; Landmark es un "
    "análisis de sensibilidad restringido a genes directamente medidos."
)

print("\nArchivos generados:")
print("- SIMILITUD_prueba_pareada_exacta.csv")
print("- SIMILITUD_prueba_pareada_detalle.csv")

print("\nFIN")