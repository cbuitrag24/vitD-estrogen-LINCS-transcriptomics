import pandas as pd

# ============================================================
# CURACION A PRIORI DE PERTURBAGENOS
# Proyecto: Estradiol vs vitamina D - LINCS L1000
# ============================================================

firmas = pd.read_csv("firmas_seleccionadas_level5.csv")
quimica = pd.read_csv("data/curation/tabla_curacion_quimica.csv")

# E2 canónico seleccionado previamente
E2_PRIMARY = ["BRD-K18910433"]

# Los registros de vitamina D se conservan por ahora para
# extracción. Su uso como análisis primario o de sensibilidad
# se definirá tras completar la verificación química.
VITD_CANDIDATES = [
    "BRD-A31158127",
    "BRD-K56429665",
    "BRD-K96390176",
    "BRD-A66229260",
    "BRD-A84327315",
    "BRD-K27316855",
]

quimica["curation_status"] = "REVIEW"

quimica.loc[
    quimica["pert_id"].isin(E2_PRIMARY),
    "curation_status"
] = "PRIMARY_E2"

quimica.loc[
    quimica["pert_id"].isin(VITD_CANDIDATES),
    "curation_status"
] = "VITD_CANDIDATE"

quimica["curation_stage"] = "PRE_EXPRESSION_ANALYSIS"

quimica.to_csv(
    "registro_curacion_preanalisis.csv",
    index=False
)

# Conservamos todas las firmas candidatas para la extracción.
# Esto evita tener que volver a acceder al GCTX de ~22 GB.
firmas.to_csv(
    "firmas_para_extraccion_level5.csv",
    index=False
)

print("============================================")
print("CURACION PREVIA AL ANALISIS DE EXPRESION")
print("============================================")

print("\nEstado de curacion:")
print(quimica["curation_status"].value_counts())

print("\nFirmas conservadas para extraccion:")
print(len(firmas))

print("\nFirmas por compuesto:")
print(firmas["pert_iname"].value_counts())

print("\nArchivos creados:")
print(" - registro_curacion_preanalisis.csv")
print(" - firmas_para_extraccion_level5.csv")

print("\nIMPORTANTE:")
print("La seleccion se realizo antes de inspeccionar")
print("los valores de expresion Level 5.")