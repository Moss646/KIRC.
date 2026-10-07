# hif_gene_sets_zscore_v1.py
# Third scoring method for the HIF panel: mean of per-gene z-scores
# (self-contained signature score). Complements ssGSEA and GSVA.
# Output: results/tables/hif_gene_sets_zscore_stats_v1.csv

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
                   "hif_gene_sets_zscore_stats_v1.csv")

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
COLLECTION = {
    "HALLMARK_HYPOXIA": "Hallmark",
    "GOBP_CELLULAR_RESPONSE_TO_DECREASED_OXYGEN_LEVELS": "GO:BP",
    "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA": "Reactome",
    "PID_HIF1_TFPATHWAY": "PID",
    "PID_HIF2PATHWAY": "PID",
    "PID_HIF1A_PATHWAY": "PID",
    "BIOCARTA_HIF_PATHWAY": "BioCarta",
    "SEMENZA_HIF1_TARGETS": "C2:CGP",
    "ELVIDGE_HYPOXIA_UP": "C2:CGP",
    "MANALO_HYPOXIA_UP": "C2:CGP",
    "WINTER_HYPOXIA_UP": "C2:CGP",
    "BUFFA_HYPOXIA_METAGENE": "C2:CGP",
    "JIANG_HYPOXIA_VIA_VHL": "C2:CGP",
}


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

    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1) + 1e-12,
                                                axis=0)

    rows = []
    scores = {}
    for s in SETS:
        genes = [g for g in all_sets[s] if g in expr.index]
        sc = z.loc[genes].mean(axis=0)
        scores[s] = sc
        u, p = mannwhitneyu(sc[s1], sc[s2])
        n1, n2 = len(s1), len(s2)
        r_rb = 2 * u / (n1 * n2) - 1
        rows.append({
            "set_name": s,
            "collection": COLLECTION[s],
            "n_genes_in_matrix": len(genes),
            "S1_mean": round(float(sc[s1].mean()), 4),
            "S2_mean": round(float(sc[s2].mean()), 4),
            "diff_S1_minus_S2": round(float(sc[s1].mean() - sc[s2].mean()), 4),
            "rank_biserial": round(float(r_rb), 4),
            "higher": "S1" if sc[s1].mean() > sc[s2].mean() else "S2",
            "p_raw": p,
        })
    df = pd.DataFrame(rows)
    df["p_bh"] = multipletests(df["p_raw"], method="fdr_bh")[1]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, index=False)

    for _, r in df.iterrows():
        sig = ("***" if r["p_bh"] < 0.001 else "**" if r["p_bh"] < 0.01
               else "*" if r["p_bh"] < 0.05 else "ns")
        print(f"  {r['set_name']:<64s} S1={r['S1_mean']:+.3f} "
              f"S2={r['S2_mean']:+.3f}  r_rb={r['rank_biserial']:+.3f}  "
              f"p_bh={r['p_bh']:.3g}  {sig}")
    print(f"wrote: {OUT}")


if __name__ == "__main__":
    main()