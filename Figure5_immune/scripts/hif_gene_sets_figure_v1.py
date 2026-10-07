# hif_gene_sets_figure_v1.py
# Summary figure for the HIF/hypoxia target gene set analysis.
# Shows rank-biserial effect size (S1 vs S2, negative = higher in S2)
# for 13 recognized gene sets across three scoring methods.
# Output: results/figures/hif_gene_sets_three_methods_v1.{png,svg}

import sys, os, io, atexit
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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

TBL = {
    "mean-z": os.path.join(REPO_ROOT, "results", "tables",
                           "hif_gene_sets_zscore_stats_v1.csv"),
    "ssGSEA": os.path.join(REPO_ROOT, "results", "tables",
                           "hif_gene_sets_stats_v1.csv"),
    "GSVA": os.path.join(REPO_ROOT, "results", "tables",
                         "hif_gene_sets_gsva_stats_v1.csv"),
}
OUT_PNG = os.path.join(REPO_ROOT, "results", "figures",
                       "hif_gene_sets_three_methods_v1.png")
OUT_SVG = os.path.join(REPO_ROOT, "results", "figures",
                       "hif_gene_sets_three_methods_v1.svg")

S2_COLOR = "#4C78A8"
S1_COLOR = "#C65A5A"
NS_COLOR = "#B9B9B9"
METHOD_MARKERS = {"mean-z": "o", "ssGSEA": "s", "GSVA": "^"}

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["svg.fonttype"] = "none"

SHORT = {
    "HALLMARK_HYPOXIA": "Hallmark Hypoxia",
    "GOBP_CELLULAR_RESPONSE_TO_DECREASED_OXYGEN_LEVELS":
        "GO:BP decr. O2 resp.",
    "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA":
        "Reactome hypoxia resp.",
    "PID_HIF1_TFPATHWAY": "PID HIF1 TF",
    "PID_HIF2PATHWAY": "PID HIF2",
    "PID_HIF1A_PATHWAY": "PID HIF1A",
    "BIOCARTA_HIF_PATHWAY": "BioCarta HIF",
    "SEMENZA_HIF1_TARGETS": "Semenza HIF1",
    "ELVIDGE_HYPOXIA_UP": "Elvidge hypoxia",
    "MANALO_HYPOXIA_UP": "Manalo hypoxia",
    "WINTER_HYPOXIA_UP": "Winter hypoxia",
    "BUFFA_HYPOXIA_METAGENE": "Buffa metagene",
    "JIANG_HYPOXIA_VIA_VHL": "Jiang via VHL",
}


def main():
    stats = {k: pd.read_csv(v) for k, v in TBL.items()}
    sets = list(stats["mean-z"]["set_name"])
    methods = ["mean-z", "ssGSEA", "GSVA"]
    offs = {"mean-z": -0.25, "ssGSEA": 0.0, "GSVA": 0.25}

    fig, ax = plt.subplots(figsize=(7.09, 4.0))
    for i, s in enumerate(sets):
        for m in methods:
            r = stats[m].set_index("set_name").loc[s]
            col = (S2_COLOR if r["rank_biserial"] < 0 else S1_COLOR)
            if r["p_bh"] >= 0.05:
                col = NS_COLOR
            ax.scatter(i + offs[m], r["rank_biserial"], s=52,
                       marker=METHOD_MARKERS[m], color=col,
                       edgecolor="black", linewidth=0.7, zorder=3)

    ax.axhline(0, color="black", lw=0.9, zorder=1)
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([SHORT[s] for s in sets], rotation=38,
                       ha="right", fontsize=8)
    ax.set_ylabel("rank-biserial (S1 vs S2)\nnegative = higher in S2",
                  fontsize=9.5)
    ax.tick_params(labelsize=8.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.7, len(sets) - 0.3)
    ax.set_ylim(-1.0, 0.9)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=S2_COLOR, markeredgecolor="black",
               label="higher in S2 (p_BH<0.05)"),
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=S1_COLOR, markeredgecolor="black",
               label="higher in S1 (p_BH<0.05)"),
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=NS_COLOR, markeredgecolor="black",
               label="ns"),
        Line2D([0], [0], marker="o", linestyle="", color="black",
               markersize=7, label="mean-z"),
        Line2D([0], [0], marker="s", linestyle="", color="black",
               markersize=7, label="ssGSEA"),
        Line2D([0], [0], marker="^", linestyle="", color="black",
               markersize=7, label="GSVA"),
    ]
    ax.legend(handles=handles, fontsize=7.5, loc="lower left",
              frameon=False, ncol=2)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=600)
    fig.savefig(OUT_SVG)
    plt.close(fig)
    print(f"wrote: {OUT_PNG}")
    print(f"wrote: {OUT_SVG}")


if __name__ == "__main__":
    main()