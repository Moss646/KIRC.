# hif_gene_sets_per_gene_check_v1.py
# Per-gene Wilcoxon test for every gene in each HIF panel set.
# Reports, per set, how many genes are significantly higher in S1 vs S2.
# Output: results/tables/hif_gene_sets_per_gene_check_v1.csv

import sys, os, io, atexit
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _hold():
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.getch()
    except Exception:
        pass


atexit.register(_hold)

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)

EXPR_PATH = os.path.join(REPO_ROOT, "data", "tcga_kirc",
                         "KIRC_expr_log2_tpm.csv")
SUB_PATH = os.path.join(REPO_ROOT, "data", "processed",
                        "subtype_assignment_balanced.csv")
GMT_H = os.path.join(REPO_ROOT, "data", "gene_sets",
                     "h.all.v2024.1.Hs.symbols.gmt")
GMT_C2 = os.path.join(REPO_ROOT, "data", "gene_sets",
                      "c2.all.v2024.1.Hs.symbols.gmt")
GMT_C5 = os.path.join(REPO_ROOT, "data", "gene_sets",
                      "c5.all.v2023.2.Hs.symbols.gmt")
OUT = os.path.join(REPO_ROOT, "results", "tables",
                   "hif_gene_sets_per_gene_check_v1.csv")

SETS = [
    "HALLMARK_HYPOXIA",
    "GOBP_CELLULAR_RESPONSE_TO_DECREASED_OXYGEN_LEVELS",
    "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
    "PID_HIF1_TFPATHWAY",
    "PID_HIF2PATHWAY",
    "PID_HIF1A_PATHWAY",
    "BIOCARTA_HIF_PATHWAY",
    "SEMENZA_HIF1_TARGETS",
    "ELVIDGE_HYPOXIA_UP",
    "MANALO_HYPOXIA_UP",
    "WINTER_HYPOXIA_UP",
    "BUFFA_HYPOXIA_METAGENE",
    "JIANG_HYPOXIA_VIA_VHL",
]


def read_gmt(path):
    sets = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def main():
    expr = pd.read_csv(EXPR_PATH, index_col=0)
    sub = pd.read_csv(SUB_PATH)
    s1 = sub.loc[sub.Subtype == "Subtype_1", "Patient"].tolist()
    s2 = sub.loc[sub.Subtype == "Subtype_2", "Patient"].tolist()

    all_sets = {}
    for p in (GMT_H, GMT_C2, GMT_C5):
        all_sets.update(read_gmt(p))

    panel = {}
    union = []
    for s in SETS:
        genes = [g for g in all_sets[s] if g in expr.index]
        panel[s] = genes
        union.extend(genes)
    union = sorted(set(union))
    print(f"union genes: {len(union)}")

    rows = []
    for g in union:
        u, p = mannwhitneyu(expr.loc[g, s1], expr.loc[g, s2])
        m1, m2 = expr.loc[g, s1].mean(), expr.loc[g, s2].mean()
        rows.append((g, p, m1 - m2))
    df = pd.DataFrame(rows, columns=["gene", "p_raw", "diff_S1_minus_S2"])
    df["p_bh"] = multipletests(df["p_raw"], method="fdr_bh")[1]
    df["direction"] = np.where(df["diff_S1_minus_S2"] > 0, "S1", "S2")
    df["sig"] = df["p_bh"] < 0.05
    gene_dir = df.set_index("gene")

    out = []
    for s in SETS:
        genes = panel[s]
        d = gene_dir.loc[genes]
        out.append({
            "set_name": s,
            "n_genes_in_matrix": len(genes),
            "n_sig": int(d["sig"].sum()),
            "n_sig_higher_S2": int(((d["sig"]) & (d["direction"] == "S2")).sum()),
            "n_sig_higher_S1": int(((d["sig"]) & (d["direction"] == "S1")).sum()),
            "frac_sig_higher_S2": round(float(
                ((d["sig"]) & (d["direction"] == "S2")).sum()
                / max(d["sig"].sum(), 1)), 3),
            "frac_nonsig": round(float((~d["sig"]).mean()), 3),
        })
    res = pd.DataFrame(out)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    res.to_csv(OUT, index=False)
    print(res.to_string(index=False))
    print(f"wrote: {OUT}")


if __name__ == "__main__":
    main()