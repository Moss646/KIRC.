"""Build msigdb_lipid_dedup_236.gmt for Figure3.

Reads the GO-deduplicated lipid gene-set table from this project's
data/gene_sets/ folder and writes the GMT file that 03_fgsea.R and
gsva_R.R read.

Input : data/gene_sets/msigdb_lipid_genesets_dedup.csv
Output: data/gene_sets/msigdb_lipid_dedup_236.gmt
"""

import os
import sys

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
GENE_SETS_DIR = os.path.join(PROJ_ROOT, "data", "gene_sets")

SRC = os.path.join(GENE_SETS_DIR, "msigdb_lipid_genesets_dedup.csv")
OUT = os.path.join(GENE_SETS_DIR, "msigdb_lipid_dedup_236.gmt")


def main():
    if not os.path.exists(SRC):
        sys.exit(f"missing input: {SRC}")

    df = pd.read_csv(SRC)
    required = {"Geneset_Name", "Description", "Genes"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"input is missing columns: {sorted(missing)}")

    os.makedirs(GENE_SETS_DIR, exist_ok=True)

    n_sets = 0
    with open(OUT, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            name = str(row["Geneset_Name"]).strip()
            desc = str(row["Description"]).strip()
            genes = [g.strip() for g in str(row["Genes"]).split("|") if g.strip()]
            if not genes:
                continue
            f.write(name + "\t" + desc + "\t" + "\t".join(genes) + "\n")
            n_sets += 1

    print(f"source: {SRC}")
    print(f"wrote : {OUT} ({n_sets} gene sets)")


if __name__ == "__main__":
    main()