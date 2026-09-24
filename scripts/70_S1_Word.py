import pandas as pd
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Cm, Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn


INPUT_FILE = Path("Supplementary_Table_S1_PUBLICATION.csv")
OUTPUT_FILE = Path("Supplementary_Table_S1_WORD.docx")

df = pd.read_csv(INPUT_FILE)

if len(df) != 12:
    raise ValueError(
        f"S1 debería contener 12 filas, pero contiene {len(df)}."
    )


def shade_cell(cell, fill="D9E2F3"):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_cell_margins(cell, value=60):
    tc_pr = cell._tc.get_or_add_tcPr()

    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)

    for margin in ["top", "start", "bottom", "end"]:
        element = tc_mar.find(qn(f"w:{margin}"))
        if element is None:
            element = OxmlElement(f"w:{margin}")
            tc_mar.append(element)

        element.set(qn("w:w"), str(value))
        element.set(qn("w:type"), "dxa")


def repeat_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def prevent_row_split(row):
    tr_pr = row._tr.get_or_add_trPr()
    cant_split = OxmlElement("w:cantSplit")
    tr_pr.append(cant_split)


def set_cell_width(cell, width_cm):
    cell.width = Cm(width_cm)

    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.first_child_found_in("w:tcW")

    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)

    tc_w.set(qn("w:w"), str(int(width_cm * 567)))
    tc_w.set(qn("w:type"), "dxa")

def format_table(table, widths, header_size=7.5, body_size=7):
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False

    for row in table.rows:
        prevent_row_split(row)

        for j, cell in enumerate(row.cells):
            set_cell_width(cell, widths[j])
            set_cell_margins(cell)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER

            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_before = Pt(0)
                paragraph.paragraph_format.space_after = Pt(0)
                paragraph.paragraph_format.line_spacing = 1

                for run in paragraph.runs:
                    run.font.name = "Arial"
                    run.font.size = Pt(body_size)

    # Header
    for cell in table.rows[0].cells:
        shade_cell(cell)

        for paragraph in cell.paragraphs:
            paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER

            for run in paragraph.runs:
                run.bold = True
                run.font.name = "Arial"
                run.font.size = Pt(header_size)

    repeat_header(table.rows[0])


def add_panel_heading(document, text):
    p = document.add_paragraph()
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(4)

    r = p.add_run(text)
    r.bold = True
    r.font.name = "Arial"
    r.font.size = Pt(9)


document = Document()

# ------------------------------------------------------------
# PAGE SETUP
# ------------------------------------------------------------

section = document.sections[0]
section.orientation = WD_ORIENT.LANDSCAPE
section.page_width, section.page_height = (
    section.page_height,
    section.page_width,
)

section.top_margin = Cm(1.4)
section.bottom_margin = Cm(1.4)
section.left_margin = Cm(1.2)
section.right_margin = Cm(1.2)

document.styles["Normal"].font.name = "Arial"
document.styles["Normal"].font.size = Pt(8)

# ------------------------------------------------------------
# TITLE
# ------------------------------------------------------------

p = document.add_paragraph()
p.paragraph_format.space_after = Pt(5)

r = p.add_run(
    "Supplementary Table S1. "
    "LINCS L1000 perturbational signatures and constituent profiles "
    "included in the primary matched analysis."
)

r.bold = True
r.font.name = "Arial"
r.font.size = Pt(10)

# ------------------------------------------------------------
# PANEL A
# ------------------------------------------------------------

add_panel_heading(
    document,
    "A. Perturbational signatures included in the primary matched analysis"
)

panel_a_columns = [
    "Cell line",
    "Experimental prefix",
    "Matched experimental block",
    "Perturbagen",
    "LINCS perturbagen ID",
    "LINCS signature ID",
    "Dose",
    "Dose unit",
    "Exposure time",
    "Time unit",
    "N constituent profiles",
]

panel_a = df[panel_a_columns]

table_a = document.add_table(
    rows=1,
    cols=len(panel_a.columns)
)

for j, column in enumerate(panel_a.columns):
    table_a.rows[0].cells[j].text = column

for _, row in panel_a.iterrows():
    cells = table_a.add_row().cells

    for j, value in enumerate(row):
        cells[j].text = "" if pd.isna(value) else str(value)

# Total usable width ≈ 27 cm
widths_a = [
    1.3,   # Cell line
    2.0,   # Experimental prefix
    2.6,   # Matched block
    2.1,   # Perturbagen
    2.8,   # LINCS pert ID
    4.4,   # LINCS signature ID
    1.0,   # Dose
    1.2,   # Dose unit
    1.5,   # Exposure time
    1.2,   # Time unit
    2.3,   # N constituent profiles
]

format_table(
    table_a,
    widths_a,
    header_size=7,
    body_size=6.8
)

# ------------------------------------------------------------
# PANEL B
# ------------------------------------------------------------

add_panel_heading(
    document,
    "B. Constituent profile identifiers contributing to each Level 5 signature"
)

# Shorter display headings ONLY.
# Source values remain unchanged.
panel_b = df[
    [
        "Cell line",
        "Perturbagen",
        "LINCS signature ID",
        "N constituent profiles",
        "Constituent profile IDs",
    ]
].copy()

panel_b.columns = [
    "Cell line",
    "Perturbagen",
    "LINCS signature ID",
    "N profiles",
    "Constituent profile IDs",
]

table_b = document.add_table(
    rows=1,
    cols=len(panel_b.columns)
)

for j, column in enumerate(panel_b.columns):
    table_b.rows[0].cells[j].text = column

for _, row in panel_b.iterrows():
    cells = table_b.add_row().cells

    for j, value in enumerate(row):
        cells[j].text = "" if pd.isna(value) else str(value)

# Give most of the page width to the long identifiers.
widths_b = [
    1.6,    # Cell line
    2.5,    # Perturbagen
    5.2,    # LINCS signature ID
    1.7,    # N profiles
    15.8,   # Constituent profile IDs
]

format_table(
    table_b,
    widths_b,
    header_size=7.5,
    body_size=7
)

# ------------------------------------------------------------
# NOTE
# ------------------------------------------------------------

p = document.add_paragraph()
p.paragraph_format.space_before = Pt(7)
p.paragraph_format.space_after = Pt(0)

r = p.add_run("Note: ")
r.bold = True
r.font.name = "Arial"
r.font.size = Pt(8)

r = p.add_run(
    "Each row represents one Level 5 MODZ perturbational signature "
    "included in the primary matched analysis. Experimental blocks "
    "were defined by cell line and LINCS experimental prefix. "
    "N constituent profiles indicates the number of constituent "
    "profiles contributing to each aggregated Level 5 signature; "
    "these profiles were not treated as independent analytical "
    "observations. LINCS perturbagen and signature identifiers are "
    "reported as provided in the GSE92742 metadata."
)

r.font.name = "Arial"
r.font.size = Pt(8)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

document.save(OUTPUT_FILE)

print()
print("FINAL REDESIGNED S1 WORD VERSION CREATED CORRECTLY")
print("------------------------------------------------")
print(f"Output file: {OUTPUT_FILE}")
print(f"Panel A rows: {len(panel_a)}")
print(f"Panel B rows: {len(panel_b)}")
print()
print("All 12 signatures and all constituent profile IDs were retained.")
print("No CSV file was modified.")