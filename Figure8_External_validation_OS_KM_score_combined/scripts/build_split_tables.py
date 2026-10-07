"""Split the multi-cohort summary CSV into three Word tables."""

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

with open(SRC, newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = [row for row in reader if row]

CONFIG = {
    "Table1_prognostic.docx": {
        "sections": ["Cohort characteristics", "Prognostic validation"],
        "title": (
            "Table 1. Prognostic validation of the 68-gene nearest-centroid (NCC) "
            "lipid-metabolic subtype classifier across TCGA-KIRC (training) and three "
            "independent cohorts."
        ),
        "notes": [
            "Abbreviations: NCC, nearest-centroid classifier; HR, hazard ratio; CI, "
            "confidence interval; C-index, Harrell's concordance index; TCGA-KIRC, "
            "The Cancer Genome Atlas kidney renal clear cell carcinoma; ICGC RECA-EU, "
            "International Cancer Genome Consortium renal cell carcinoma European cohort; "
            "CPTAC, Clinical Proteomic Tumor Analysis Consortium.",
            "S1 and S2 denote the two lipid-metabolic subtypes "
            "(S1 = lipid-catabolic / FAO-low; S2 = lipid-anabolic / FAO-high).",
            "Prognostic (survival) validation was not performed in the CPTAC cohort (-).",
            "Stage-stratified multivariable Cox models were adjusted for pathological "
            "stage, grade, and age; complete-case analysis required non-missing stage, "
            "grade, and age.",
            "-, not assessed / not available.",
        ],
    },
    "TableS1_classifier.docx": {
        "sections": ["Classifier robustness"],
        "title": (
            "Supplementary Table S1. Transferability of the 68-gene NCC classifier "
            "across external cohorts, assessed by Pearson correlation of subtype scores "
            "versus TCGA-KIRC."
        ),
        "notes": [
            "Abbreviations: NCC, nearest-centroid classifier; TCGA-KIRC, The Cancer "
            "Genome Atlas kidney renal clear cell carcinoma; ICGC RECA-EU, International "
            "Cancer Genome Consortium renal cell carcinoma European cohort; CPTAC, "
            "Clinical Proteomic Tumor Analysis Consortium.",
            "Pearson correlation and 95% CI were computed between per-sample NCC subtype "
            "scores in each external cohort and the corresponding TCGA-KIRC scores.",
            "Missing genes were those absent from each cohort's expression matrix and "
            "therefore excluded from score computation.",
            "-, not assessed / not available.",
        ],
    },
    "TableS2_biological.docx": {
        "sections": ["Biological validation"],
        "title": (
            "Supplementary Table S2. Cross-cohort biological validation of the 14-gene "
            "lipid-metabolic panel, assessed by Mann-Whitney U test with "
            "Benjamini-Hochberg correction (FDR < 0.05) in each external cohort."
        ),
        "notes": [
            "Direction concordance = number of the 14 lipid-metabolic genes whose "
            "S1-vs-S2 expression direction matched the TCGA-KIRC direction. "
            "BH-sig genes = genes with FDR < 0.05 by Mann-Whitney U test. "
            "Pearson r and 95% CI compare per-gene effect sizes versus TCGA-KIRC.",
            "Abbreviations: MWU, Mann-Whitney U test; BH, Benjamini-Hochberg; "
            "TCGA-KIRC, The Cancer Genome Atlas kidney renal clear cell carcinoma; "
            "ICGC RECA-EU, International Cancer Genome Consortium renal cell carcinoma "
            "European cohort; CPTAC, Clinical Proteomic Tumor Analysis Consortium.",
            "-, not assessed / not available.",
        ],
    },
}


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


def build(rows_subset, title, notes, outname):
    doc = Document()
    doc.styles["Normal"].font.name = "Arial"
    doc.styles["Normal"].font.size = Pt(9)

    t = doc.add_paragraph()
    t.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = t.add_run(title)
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

    groups = []
    for row in rows_subset:
        sec = row[0]
        if not groups or groups[-1][0] != sec:
            groups.append([sec, []])
        groups[-1][1].append(row)

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

    for n in notes:
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        rn = p.add_run(n)
        rn.font.name = "Arial"
        rn.font.size = Pt(8)

    out = os.path.join(TABLES_DIR, outname)
    doc.save(out)
    return out


for outname, cfg in CONFIG.items():
    subset = [row for row in rows if row[0] in cfg["sections"]]
    out = build(subset, cfg["title"], cfg["notes"], outname)
    print(f"saved: {out}  rows: {len(subset)}  sections: {cfg['sections']}")