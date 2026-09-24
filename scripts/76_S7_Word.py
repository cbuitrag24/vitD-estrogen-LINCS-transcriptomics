import pandas as pd
from docx import Document
from docx.shared import Cm, Pt
from docx.enum.section import WD_ORIENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

FILE_A = "Supplementary_Table_S7A.csv"
FILE_B = "Supplementary_Table_S7B.csv"
OUTPUT_FILE = "Supplementary_Table_S7_WORD.docx"

TITLE = (
    "Supplementary Table S7. Context robustness of the eight core convergent "
    "cell-cycle genes across estradiol, calcitriol and calcipotriol."
)

NOTE = (
    "Note: Table S7A summarises the behaviour of the eight robust convergent "
    "genes within each treatment across the four matched experimental contexts, "
    "including the numbers of negative, positive and neutral contexts together "
    "with mean, median and range of MODZ values. Table S7B reports the integrated "
    "cross-treatment summary, including the mean direction of regulation, "
    "treatment-specific mean MODZ values, context-consistency classification and "
    "whether each gene satisfied the robustness criterion simultaneously across "
    "all three treatments."
)

# ------------------------------------------------------------
# Read audited source files
# ------------------------------------------------------------

a = pd.read_csv(FILE_A)
b = pd.read_csv(FILE_B)

if len(a) != 39:
    raise ValueError(f"Unexpected S7A row count: {len(a)}. Expected 39.")

if len(b) != 13:
    raise ValueError(f"Unexpected S7B row count: {len(b)}. Expected 13.")

# ------------------------------------------------------------
# Publication copies (editorial headings only)
# ------------------------------------------------------------

a_pub = a.rename(columns={
    "Treatment": "Treatment",
    "Gene": "Gene",
    "N_contexts": "Contexts",
    "N_negative": "Negative",
    "N_positive": "Positive",
    "N_zero": "Neutral",
    "Mean_MODZ": "Mean MODZ",
    "Median_MODZ": "Median MODZ",
    "Min_MODZ": "Minimum MODZ",
    "Max_MODZ": "Maximum MODZ",
    "Context_consistency": "Context consistency",
}).copy()

b_pub = b.rename(columns={
    "Gene": "Gene",
    "Mean_direction_across_treatments": "Mean direction",
    "Estradiol_mean_MODZ": "Estradiol",
    "Calcitriol_mean_MODZ": "Calcitriol",
    "Calcipotriol_mean_MODZ": "Calcipotriol",
    "Estradiol_context_consistency": "Estradiol consistency",
    "Calcitriol_context_consistency": "Calcitriol consistency",
    "Calcipotriol_context_consistency": "Calcipotriol consistency",
    "Context_robust_all_3_treatments": "Robust in all treatments",
}).copy()

# Editorial translation of categorical labels in the Word copy only
translation_map = {
    "predominio negativo 3de4": "Predominantly negative (3/4)",
    "negativo 4de4": "Negative (4/4)",
    "heterogéneo": "Heterogeneous",
    "heterogeneo": "Heterogeneous",
    "negativa en los 3": "Negative in all 3",
}

a_pub["Context consistency"] = (
    a_pub["Context consistency"]
    .replace(translation_map)
)

for column in [
    "Mean direction",
    "Estradiol consistency",
    "Calcitriol consistency",
    "Calcipotriol consistency",
]:
    b_pub[column] = b_pub[column].replace(translation_map)
b_pub["Robust in all treatments"] = (
    b_pub["Robust in all treatments"]
    .replace({True: "Yes", False: "No"})
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def display_value(value):
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)

def shade(cell):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), "D9E2F3")
    tc_pr.append(shd)

def repeat_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    hdr = OxmlElement("w:tblHeader")
    hdr.set(qn("w:val"), "true")
    tr_pr.append(hdr)

def add_heading(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run(text)
    r.bold = True
    r.font.name = "Arial"
    r.font.size = Pt(9)

def add_table(doc, df):
    table = doc.add_table(rows=1, cols=len(df.columns))
    table.style = "Table Grid"
    table.autofit = True

    repeat_header(table.rows[0])

    for j, col in enumerate(df.columns):
        c = table.rows[0].cells[j]
        c.text = str(col)
        c.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        shade(c)
        for p in c.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for r in p.runs:
                r.bold = True
                r.font.name = "Arial"
                r.font.size = Pt(7)

    for _, row in df.iterrows():
        cells = table.add_row().cells
        for j, val in enumerate(row):
            cells[j].text = display_value(val)
            cells[j].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cells[j].paragraphs:
                for r in p.runs:
                    r.font.name = "Arial"
                    r.font.size = Pt(7)

# ------------------------------------------------------------
# Word document
# ------------------------------------------------------------

doc = Document()

sec = doc.sections[0]
sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width = Cm(29.7)
sec.page_height = Cm(21.0)
sec.top_margin = Cm(1.3)
sec.bottom_margin = Cm(1.3)
sec.left_margin = Cm(1.3)
sec.right_margin = Cm(1.3)

doc.styles["Normal"].font.name = "Arial"
doc.styles["Normal"].font.size = Pt(8)

p = doc.add_paragraph()
r = p.add_run(TITLE)
r.bold = True
r.font.name = "Arial"
r.font.size = Pt(10)

add_heading(doc, "S7A. Context-specific robustness within each treatment")
add_table(doc, a_pub)

doc.add_page_break()

add_heading(doc, "S7B. Integrated cross-treatment robustness summary")
add_table(doc, b_pub)

p = doc.add_paragraph()
p.paragraph_format.space_before = Pt(6)
r = p.add_run(NOTE)
r.font.name = "Arial"
r.font.size = Pt(8)

doc.save(OUTPUT_FILE)

print("S7 WORD VERSION CREATED CORRECTLY")
print("---------------------------------")
print(f"Output file: {OUTPUT_FILE}")
print(f"S7A rows: {len(a_pub)}")
print(f"S7B rows: {len(b_pub)}")
print()
print("The audited S7A and S7B CSV files were NOT modified.")