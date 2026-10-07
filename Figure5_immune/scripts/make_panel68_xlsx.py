"""Build Supplementary_Table_NCC_68_genes.xlsx for Figure5.

Reads the 68 gene names from the NCC classifier JSON and the Cox results
from the strict selection CSV, adds a Role column derived from HR, and
writes a formatted Excel file with the columns that
panel68_cell_source_binning_v1.py expects:

    data/processed/Supplementary_Table_NCC_68_genes.xlsx

Columns (A-G):
    Gene, HR, log2HR, p_value, FDR, Role, Functional_Layer

The reading script uses column A (Gene), F (Role) and G (Functional_Layer),
starting from row 4.
"""

import json
import os
import sys

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA = os.path.join(PROJ_ROOT, "data")
PROCESSED = os.path.join(DATA, "processed")

SEL_CSV = os.path.join(PROCESSED, "selected_genes_strict.csv")
CLF_JSON = os.path.join(PROCESSED, "subtype_classifier.json")
OUT = os.path.join(PROCESSED, "Supplementary_Table_NCC_68_genes.xlsx")


def write_cell(cell, col, v):
    if pd.isna(v):
        cell.value = None
        return
    if col == "HR":
        cell.value = float(v)
        cell.number_format = "0.0000"
    elif col == "log2HR":
        cell.value = float(v)
        cell.number_format = "0.0000"
    elif col == "p_value":
        cell.value = float(v)
        cell.number_format = "0.00E+00"
    elif col == "FDR":
        cell.value = float(v)
        cell.number_format = "0.00E+00"
    else:
        cell.value = v


def main():
    if not os.path.exists(SEL_CSV):
        sys.exit(f"missing input: {SEL_CSV}")
    if not os.path.exists(CLF_JSON):
        sys.exit(f"missing input: {CLF_JSON}")

    with open(CLF_JSON) as f:
        clf = json.load(f)
    genes68 = list(clf["genes"])
    if len(genes68) != 68:
        print(f"WARNING: classifier JSON has {len(genes68)} genes, expected 68")

    df = pd.read_csv(SEL_CSV)
    if "gene" not in df.columns or "HR" not in df.columns:
        sys.exit("selected_genes_strict.csv must contain 'gene' and 'HR'")

    missing = [g for g in genes68 if g not in set(df["gene"])]
    if missing:
        print(f"WARNING: {len(missing)} genes not in selected_genes_strict.csv:")
        for g in missing:
            print(f"  {g}")

    df = df.set_index("gene").reindex(genes68).reset_index()
    df = df.rename(columns={"gene": "Gene"})

    df["Role"] = df["HR"].apply(
        lambda hr: "Risk" if pd.notna(hr) and hr > 1 else "Protective"
    )

    out = df[["Gene", "HR", "log2HR", "p_value", "FDR",
              "Role", "functional_layer"]].copy()
    out = out.rename(columns={"functional_layer": "Functional_Layer"})

    headers = ["Gene", "HR", "log2HR", "p_value", "FDR",
               "Role", "Functional_Layer"]

    wb = Workbook()
    ws = wb.active
    ws.title = "NCC Classifier"

    ws.merge_cells("A1:G1")
    ws["A1"] = ("Supplementary Table. 68 lipid metabolism-related genes "
                "selected for NCC classifier construction.")
    ws["A1"].font = Font(bold=True, size=12)
    ws["A1"].alignment = Alignment(wrap_text=True)
    ws.row_dimensions[1].height = 25

    descs = {
        "A": "Gene: Gene symbol (HGNC)",
        "B": "HR: Hazard ratio from univariate Cox regression",
        "C": "log2HR: log2(Hazard ratio)",
        "D": "p_value: Nominal P-value from univariate Cox regression",
        "E": "FDR: Benjamini-Hochberg false discovery rate",
        "F": "Role: Risk (HR > 1) or Protective (HR < 1)",
        "G": "Functional_Layer: Major lipid metabolism functional category "
             "assigned by MSigDB gene set membership",
    }
    for col, desc in descs.items():
        ws[f"{col}2"] = desc
        ws[f"{col}2"].font = Font(italic=True, size=9, color="555555")
        ws[f"{col}2"].alignment = Alignment(wrap_text=True)
    ws.row_dimensions[2].height = 40

    for i, h in enumerate(headers, 1):
        c = ws.cell(row=3, column=i, value=h)
        c.font = Font(bold=True, size=10)
        c.alignment = Alignment(horizontal="center")

    for r_idx, (_, row) in enumerate(out.iterrows(), 4):
        for c_idx, col in enumerate(headers, 1):
            write_cell(ws.cell(row=r_idx, column=c_idx), col, row[col])

    ws.column_dimensions["A"].width = 14
    for c in "BCDEFG":
        ws.column_dimensions[c].width = 22

    wb.save(OUT)
    print(f"wrote {OUT} ({len(out)} genes)")


if __name__ == "__main__":
    main()