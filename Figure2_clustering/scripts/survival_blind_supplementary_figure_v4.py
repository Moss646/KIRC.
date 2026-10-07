# -*- coding: utf-8 -*-
"""
survival_blind_supplementary_figure_v4.py
=========================================
v4: reads the merged analysis outputs (survival_blind_analysis_v1_*) from
survival_blind_analysis_v1.py instead of the two legacy upstream scripts.
Layout, publication size (180 mm) and the Panel C "Consensus index" colorbar
are identical to v3; v1-v3 scripts and figures untouched.

Panels
------
A. PAC (all pairs, Monti original definition) across k = 2-6.
B. Mean silhouette width in expression space (Euclidean, z-scored lipid
   genes) across k = 2-6.
C. Survival-blind consensus matrix at k = 2, ordered by blind cluster.
D. Sankey diagram: survival-blind clusters vs published S1/S2 subtypes.

Outputs
-------
- results/figures/survival_blind_supplementary_figure_v4.png / .svg
"""

import io
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as mpatches
from scipy.stats import chi2_contingency
from sklearn.metrics import adjusted_rand_score, cohen_kappa_score

warnings.filterwarnings("ignore")
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["xtick.major.width"] = 0.8
plt.rcParams["ytick.major.width"] = 0.8

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(HERE)

EXPR_FILE = os.path.join(PROJ_ROOT, "data", "KIRC_expr_log2_tpm.csv")
SUB_FILE = os.path.join(PROJ_ROOT, "data",
                        "subtype_assignment_balanced.csv")
OUT_DIR = os.path.join(PROJ_ROOT, "results")
FIG_DIR = os.path.join(OUT_DIR, "figures")

STATS_FILE = os.path.join(OUT_DIR, "survival_blind_analysis_v1_stats.csv")
CONS_K2 = os.path.join(OUT_DIR,
                       "survival_blind_analysis_v1_consensus_k2.npy")
ASSIGN_FILE = os.path.join(OUT_DIR,
                           "survival_blind_analysis_v1_assignment.csv")

S1_COLOR = "#C65A5A"
S2_COLOR = "#4C78A8"
GRAY = "#B8B8B8"
NODE_GRAY = "#5A5A5A"
HIGHLIGHT = "#C65A5A"
DEFAULT_BAR = "#9AA5B1"

FIG_W_IN = 7.087
FIG_H_IN = 6.10
TITLE_FS = 10
LABEL_FS = 9
TICK_FS = 8
BARLABEL_FS = 7.5
SANKEY_FS = 8


def sankey_ribbon(ax, x0, y0_top, y0_bot, x1, y1_top, y1_bot, color, alpha):
    xm = 0.5 * (x0 + x1)
    verts = [(x0, y0_top), (xm, y0_top), (xm, y1_top), (x1, y1_top),
             (x1, y1_bot), (xm, y1_bot), (xm, y0_bot), (x0, y0_bot),
             (x0, y0_top)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), facecolor=color,
                                    edgecolor="none", alpha=alpha, zorder=1))


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    st = pd.read_csv(STATS_FILE).sort_values("k")
    ks = st["k"].to_numpy()
    pac = st["PAC_allpairs"].to_numpy()
    sil = st["mean_silhouette_expression"].to_numpy()

    assign = pd.read_csv(ASSIGN_FILE)
    expr_cols = pd.read_csv(EXPR_FILE, index_col=0, nrows=0).columns.tolist()
    sub = pd.read_csv(SUB_FILE).set_index("Patient")["Subtype"]
    sample_order = [s for s in expr_cols if s in sub.index]
    a = assign.set_index("Patient").loc[sample_order]
    y_sub = a["Published"].to_numpy()
    blind = a["Survival_blind_cluster"].to_numpy()
    assert set(np.unique(blind)) == {"S1-blind", "S2-blind"}
    m1 = blind == "S1-blind"
    m2 = blind == "S2-blind"
    n_samp = len(sample_order)
    print(f"samples: {n_samp} (blind1={m1.sum()}, blind2={m2.sum()}, "
          f"S1={(y_sub == 'Subtype_1').sum()}, "
          f"S2={(y_sub == 'Subtype_2').sum()})")

    blind_mapped = np.where(m1, "Subtype_1", "Subtype_2")
    ari = adjusted_rand_score(y_sub, blind_mapped)
    kappa = cohen_kappa_score(y_sub, blind_mapped)
    ct = pd.crosstab(pd.Series(y_sub, name="Published"),
                     pd.Series(blind, name="Blind"))
    ct = ct[["S1-blind", "S2-blind"]]
    n11 = int(ct.loc["Subtype_1", "S1-blind"])
    n12 = int(ct.loc["Subtype_1", "S2-blind"])
    n21 = int(ct.loc["Subtype_2", "S1-blind"])
    n22 = int(ct.loc["Subtype_2", "S2-blind"])
    n_conc = n11 + n22
    conc = n_conc / n_samp
    print(ct.to_string())
    print(f"concordance={conc:.4f} ({n_conc}/{n_samp}) ARI={ari:.4f} "
          f"kappa={kappa:.4f}")

    C = np.load(CONS_K2)
    assert C.shape == (n_samp, n_samp), \
        f"consensus {C.shape} vs samples {n_samp}"
    order = np.argsort(~m1, kind="stable")

    fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.15],
                          hspace=0.48, wspace=0.34)

    axA = fig.add_subplot(gs[0, 0])
    colors_a = [HIGHLIGHT if k == 2 else DEFAULT_BAR for k in ks]
    axA.bar(ks, pac, width=0.55, color=colors_a, edgecolor="black", lw=0.6)
    for k, v in zip(ks, pac):
        axA.text(k, v + 0.012, f"{v:.2f}",
                 ha="center", va="bottom", fontsize=BARLABEL_FS)
    axA.set_xticks(ks)
    axA.set_xlabel("Number of clusters (k)", fontsize=LABEL_FS)
    axA.set_ylabel("PAC (all pairs)", fontsize=LABEL_FS)
    axA.set_ylim(0, 0.55)
    axA.tick_params(labelsize=TICK_FS)
    axA.set_title("A", fontsize=TITLE_FS, loc="left", fontweight="bold")

    axB = fig.add_subplot(gs[0, 1])
    colors_b = [HIGHLIGHT if k == 2 else DEFAULT_BAR for k in ks]
    axB.bar(ks, sil, width=0.55, color=colors_b, edgecolor="black", lw=0.6)
    for k, v in zip(ks, sil):
        axB.text(k, v + 0.006, f"{v:.2f}", ha="center", va="bottom",
                 fontsize=BARLABEL_FS)
    axB.set_xticks(ks)
    axB.set_xlabel("Number of clusters (k)", fontsize=LABEL_FS)
    axB.set_ylabel("Mean silhouette (expression)", fontsize=LABEL_FS)
    axB.set_ylim(0, 0.32)
    axB.tick_params(labelsize=TICK_FS)
    axB.set_title("B", fontsize=TITLE_FS, loc="left", fontweight="bold")

    axC = fig.add_subplot(gs[1, 0])
    im = axC.imshow(C[np.ix_(order, order)], cmap="Blues", vmin=0, vmax=1,
                    interpolation="nearest", aspect="equal")
    n1 = int(m1.sum())
    axC.axhline(n1 - 0.5, color="black", lw=0.8)
    axC.axvline(n1 - 0.5, color="black", lw=0.8)
    axC.set_xticks([])
    axC.set_yticks([])
    axC.set_xlabel("Samples (ordered by blind cluster)", fontsize=LABEL_FS)
    axC.set_ylabel("Samples", fontsize=LABEL_FS)
    axC.set_title("C", fontsize=TITLE_FS, loc="left", fontweight="bold")
    cbar = fig.colorbar(im, ax=axC, fraction=0.05, pad=0.04,
                        ticks=[0, 0.5, 1])
    cbar.ax.tick_params(labelsize=TICK_FS, width=0.8)
    cbar.outline.set_linewidth(0.8)
    cbar.set_label("Consensus index", fontsize=LABEL_FS)

    axD = fig.add_subplot(gs[1, 1])
    axD.set_xlim(0, 1)
    axD.set_ylim(-0.24, 1.18)
    axD.axis("off")
    xl, xr, w = 0.17, 0.76, 0.08

    left_nodes = [("Blind cluster 1", m1.sum(), NODE_GRAY),
                  ("Blind cluster 2", m2.sum(), NODE_GRAY)]
    right_nodes = [("S1", (y_sub == "Subtype_1").sum(), S1_COLOR),
                   ("S2", (y_sub == "Subtype_2").sum(), S2_COLOR)]

    def spans(counts_seq):
        edges = np.cumsum([0] + list(counts_seq)) / n_samp
        return [(1 - edges[i + 1], 1 - edges[i])
                for i in range(len(counts_seq))]

    left_spans = spans([m1.sum(), m2.sum()])
    right_spans = spans([(y_sub == "Subtype_1").sum(),
                         (y_sub == "Subtype_2").sum()])
    left_top = [0, int(m1.sum())]
    right_top = [0, int((y_sub == "Subtype_1").sum())]
    flow_draw = [(0, 0, n11, S1_COLOR), (1, 0, n21, GRAY),
                 (0, 1, n12, GRAY), (1, 1, n22, S2_COLOR)]
    left_used = [0, 0]
    right_used = [0, 0]
    label_x = {(1, 0): 0.40, (0, 1): 0.60}
    for li, ri, cntv, color in flow_draw:
        lt = 1 - (left_top[li] + left_used[li] + cntv) / n_samp
        lb = 1 - (left_top[li] + left_used[li]) / n_samp
        rt = 1 - (right_top[ri] + right_used[ri] + cntv) / n_samp
        rb = 1 - (right_top[ri] + right_used[ri]) / n_samp
        sankey_ribbon(axD, xl + w, lt, lb, xr, rt, rb, color, alpha=0.65)
        axD.text(label_x.get((li, ri), 0.5),
                 (lt + lb + rt + rb) / 4, str(cntv),
                 ha="center", va="center", fontsize=SANKEY_FS, zorder=3,
                 color="#333333",
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                           alpha=0.8))
        left_used[li] += cntv
        right_used[ri] += cntv
    for (name, cntv, color), (t, b) in zip(left_nodes, left_spans):
        axD.add_patch(mpatches.Rectangle(
            (xl, b), w, t - b, facecolor=color,
            edgecolor="black", lw=0.6, zorder=2))
        axD.text(xl - 0.03, (t + b) / 2, f"{name}\n(n = {int(cntv)})",
                 ha="right", va="center", fontsize=SANKEY_FS)
    for (name, cntv, color), (t, b) in zip(right_nodes, right_spans):
        axD.add_patch(mpatches.Rectangle(
            (xr, b), w, t - b, facecolor=color,
            edgecolor="black", lw=0.6, zorder=2))
        axD.text(xr + w + 0.03, (t + b) / 2, f"{name}\n(n = {int(cntv)})",
                 ha="left", va="center", fontsize=SANKEY_FS)
    axD.text(xl + w / 2, 1.06, "Survival-blind\nclusters", ha="center",
             va="bottom", fontsize=SANKEY_FS)
    axD.text(xr + w / 2, 1.06, "Published\nsubtypes", ha="center",
             va="bottom", fontsize=SANKEY_FS)
    axD.text(0.5, -0.12,
             f"Concordance = {conc * 100:.1f}% ({n_conc}/{n_samp})\n"
             f"ARI = {ari:.2f}, Cohen's $\\kappa$ = {kappa:.2f}",
             ha="center", va="top", fontsize=SANKEY_FS,
             transform=axD.transAxes)
    axD.set_title("D", fontsize=TITLE_FS, loc="left",
                  fontweight="bold", pad=24)

    fig.savefig(os.path.join(
        FIG_DIR, "survival_blind_supplementary_figure_v4.png"),
        dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(
        FIG_DIR, "survival_blind_supplementary_figure_v4.svg"),
        bbox_inches="tight")
    plt.close(fig)
    print("figure written")
    return 0


if __name__ == "__main__":
    sys.exit(main())