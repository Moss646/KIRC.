"""Build the merged main Word table from the multi-cohort summary CSV."""

import csv
import os

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

TABLES_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "results", "tables",
)
SRC = os.path.join(TABLES_DIR, "NCC_multi_cohort_validation_summary.csv")
OUT = os.path.join(TABLES_DIR, "NCC_validation_Table1.docx")

with open(SRC, newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = [row for row in reader if row]

groups = []
for row in rows:
    sec = row[0]
    if not groups or groups[-1][0] != sec:
        groups.append([sec, []])
    groups[-1][1].append(row)


def set_cell(cell, text):
    cell.text = ""
    text = text.strip()
    if "e" in text or "E" in text:
        try:
            sep = "e" if "e" in text else "E"
            mantissa, exp = text.split(sep)
            mantissa = mantissa.strip()
            exp = int(exp.strip())
            p = cell.paragraphs[0]
            p.add_run(f"{mantissa} \u00d7 10")
            run = p.add_run(str(exp))
            run.font.superscript = True
            return
        except ValueError:
            pass
    cell.text = text


def style_cell_text(cell, bold=False, size=9, align=None):
    for p in cell.paragraphs:
        if align is not None:
            p.alignment = align
        for rn in p.runs:
            rn.font.name = "Arial"
            rn.font.size = Pt(size)
            rn.font.bold = bold


doc = Document()
normal = doc.styles["Normal"]
normal.font.name = "Arial"
normal.font.size = Pt(9)

title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.LEFT
run = title.add_run(
    "Table 1. Multi-cohort external validation of the 68-gene nearest-centroid "
    "(NCC) lipid-metabolic subtype classifier across TCGA-KIRC (training) and three "
    "independent cohorts."
)
run.bold = True
run.font.name = "Arial"
run.font.size = Pt(10)

n_cohorts = len(header) - 2
table = doc.add_table(rows=1, cols=2 + n_cohorts)
table.style = "Table Grid"
table.alignment = WD_TABLE_ALIGNMENT.CENTER

hdr = table.rows[0].cells
hdr[0].text = "Category"
hdr[1].text = "Metric"
for i, c in enumerate(header[2:]):
    hdr[2 + i].text = c

for c in hdr:
    style_cell_text(c, bold=True, size=9, align=WD_ALIGN_PARAGRAPH.CENTER)
    tcPr = c._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:fill"), "D9E2F3")
    tcPr.append(shd)

body_idx = 0
for sec_name, sec_rows in groups:
    start = body_idx
    for row in sec_rows:
        cells = table.add_row().cells
        set_cell(cells[1], row[1])
        for i, val in enumerate(row[2:]):
            set_cell(cells[2 + i], val)
            style_cell_text(cells[2 + i], size=9, align=WD_ALIGN_PARAGRAPH.CENTER)
        style_cell_text(cells[1], size=9, align=WD_ALIGN_PARAGRAPH.LEFT)
        body_idx += 1
    end = body_idx - 1
    top = start + 1
    bot = end + 1
    if end > start:
        merged = table.cell(top, 0).merge(table.cell(bot, 0))
        merged.text = sec_name
        style_cell_text(merged, bold=True, size=9, align=WD_ALIGN_PARAGRAPH.LEFT)
    else:
        c = table.cell(top, 0)
        c.text = sec_name
        style_cell_text(c, bold=True, size=9, align=WD_ALIGN_PARAGRAPH.LEFT)

notes = [
    "Abbreviations: NCC, nearest-centroid classifier; HR, hazard ratio; CI, "
    "confidence interval; C-index, Harrell's concordance index; MWU, Mann-Whitney U "
    "test; BH, Benjamini-Hochberg; TCGA-KIRC, The Cancer Genome Atlas kidney renal "
    "clear cell carcinoma; ICGC RECA-EU, International Cancer Genome Consortium "
    "renal cell carcinoma European cohort; CPTAC, Clinical Proteomic Tumor Analysis "
    "Consortium.",
    "S1 and S2 denote the two lipid-metabolic subtypes "
    "(S1 = lipid-catabolic / FAO-low; S2 = lipid-anabolic / FAO-high).",
    "Prognostic (survival) validation was not performed in the CPTAC cohort (-).",
    "Stage-stratified multivariable Cox models were adjusted for pathological stage, "
    "grade, and age; complete-case analysis required non-missing stage, grade, and age.",
    "Direction concordance denotes the number of genes (out of all available genes of "
    "the 14-gene lipid panel in each cohort; protein abundance for CPTAC) whose "
    "direction agrees with the TCGA-KIRC reference direction. BH-significant gene "
    "counts denote the number of these available genes showing significant "
    "differential expression by Mann-Whitney U test with Benjamini-Hochberg "
    "correction (FDR < 0.05).",
    "-, not assessed / not available.",
]
for n in notes:
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    rn = p.add_run(n)
    rn.font.name = "Arial"
    rn.font.size = Pt(8)

doc.save(OUT)
print(f"saved: {OUT}")
print(f"rows: {len(rows)}  groups: {[g[0] for g in groups]}")