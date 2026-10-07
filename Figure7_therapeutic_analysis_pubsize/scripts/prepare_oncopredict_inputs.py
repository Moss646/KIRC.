"""Create kirc_expr_matched.csv and kirc_subtypes.csv for the Figure7 R script."""

import os
import sys

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJ_ROOT, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")

CANDIDATE_DIRS = [DATA_DIR, RAW_DIR]


def find_input(name):
    for d in CANDIDATE_DIRS:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


OUT_EXPR = os.path.join(RAW_DIR, "kirc_expr_matched.csv")
OUT_SUB = os.path.join(RAW_DIR, "kirc_subtypes.csv")


def main():
    expr_src = find_input("KIRC_expr_log2_tpm.csv")
    sub_src = find_input("subtype_assignment_balanced.csv")

    if expr_src is None:
        sys.exit("missing input: KIRC_expr_log2_tpm.csv (looked in data/ and data/raw/)")
    if sub_src is None:
        sys.exit("missing input: subtype_assignment_balanced.csv (looked in data/ and data/raw/)")

    os.makedirs(RAW_DIR, exist_ok=True)

    expr = pd.read_csv(expr_src, index_col=0)
    expr.to_csv(OUT_EXPR)
    print(f"wrote {OUT_EXPR} ({expr.shape[0]} genes x {expr.shape[1]} samples)")

    sub = pd.read_csv(sub_src)
    sub.columns = ["X", "Subtype"]
    sub.to_csv(OUT_SUB, index=False)
    print(f"wrote {OUT_SUB} ({len(sub)} samples)")


if __name__ == "__main__":
    main()