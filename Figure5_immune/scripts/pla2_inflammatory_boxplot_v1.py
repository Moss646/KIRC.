# pla2_inflammatory_boxplot_v1.py
# Immune/inflammatory PLA2-related gene set: S1 vs S2 expression + set score.
# 4 genes: PLA2G2A, PLA2G4F, PLA2G1B, PLA2G2D.
# Output: results/tables/pla2_inflammatory_stats_v1.csv
#         results/figures/pla2_inflammatory_S1S2_v1.{png,svg}

import sys, os, io, atexit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_ROOT = os.environ.get('MEL_DATA_ROOT', PROJ_ROOT)

EXPR_CSV = os.path.join(DATA_ROOT, "data", "tcga_kirc",
                        "KIRC_expr_log2_tpm.csv")
SUB_CSV = os.path.join(DATA_ROOT, "data", "processed",
                       "subtype_assignment_balanced.csv")
TAB_OUT = os.path.join(PROJ_ROOT, "results", "tables",
                       "pla2_inflammatory_stats_v1.csv")
FIG_OUT = os.path.join(PROJ_ROOT, "results", "figures",
                       "pla2_inflammatory_S1S2_v1")

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

GENES = ["PLA2G2A", "PLA2G4F", "PLA2G1B", "PLA2G2D"]
DISPLAY = {g: g for g in GENES}

expr = pd.read_csv(EXPR_CSV, index_col=0)
sub = pd.read_csv(SUB_CSV)
sm = sub.set_index("Patient")["Subtype"]
common = [c for c in expr.columns if c in sm.index]
s1_ids = [c for c in common if sm[c] == "Subtype_1"]
s2_ids = [c for c in common if sm[c] == "Subtype_2"]

missing = [g for g in GENES if g not in expr.index]
if missing:
    raise SystemExit(f"genes missing from matrix: {missing}")

rows, pvals = [], []
for g in GENES:
    v1 = expr.loc[g, s1_ids].astype(float)
    v2 = expr.loc[g, s2_ids].astype(float)
    p = mannwhitneyu(v2, v1, alternative="two-sided").pvalue
    rows.append({"gene": g, "S1_mean": v1.mean(), "S2_mean": v2.mean(),
                 "diff_S2_minus_S1": v2.mean() - v1.mean(), "MWU_p": p})
    pvals.append(p)

# Set score: mean of per-gene z-scores computed across all samples.
z = expr.loc[GENES, common].astype(float).T
z = (z - z.mean()) / z.std(ddof=0)
score = z.mean(axis=1)
s1_sc, s2_sc = score.loc[s1_ids], score.loc[s2_ids]
p_score = mannwhitneyu(s2_sc, s1_sc, alternative="two-sided").pvalue

rows.append({"gene": "SET_SCORE_mean_z", "S1_mean": s1_sc.mean(),
             "S2_mean": s2_sc.mean(),
             "diff_S2_minus_S1": s2_sc.mean() - s1_sc.mean(),
             "MWU_p": p_score})
pvals.append(p_score)

res = pd.DataFrame(rows)
res["BH_FDR"] = multipletests(np.array(pvals), method="fdr_bh")[1]
os.makedirs(os.path.dirname(TAB_OUT), exist_ok=True)
res.to_csv(TAB_OUT, index=False)

C1, C2 = "#C65A5A", "#4C78A8"
fig, (ax, ax2) = plt.subplots(
    1, 2, figsize=(7.09, 3.4),
    gridspec_kw={"width_ratios": [4, 1.6], "wspace": 0.30})
rng = np.random.default_rng(1)


def paired_boxes(ax, values, labels, ylabel):
    for k, (v1, v2, fdr) in enumerate(values):
        for off, vals, col in ((-0.18, v1, C1), (0.18, v2, C2)):
            ax.boxplot([vals], positions=[k + off], widths=0.30,
                       patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", lw=1.0),
                       boxprops=dict(facecolor=col, alpha=0.55, lw=0.8),
                       whiskerprops=dict(lw=0.8),
                       capprops=dict(lw=0.8))
            jitter = rng.normal(k + off, 0.05, size=len(vals))
            ax.scatter(jitter, vals, s=2.5, color=col, alpha=0.45, lw=0)
        if fdr < 0.001:
            st = "***"
        elif fdr < 0.01:
            st = "**"
        elif fdr < 0.05:
            st = "*"
        else:
            st = "ns"
        ymax = max(v1.max(), v2.max())
        ymin = min(v1.min(), v2.min())
        span = ymax - ymin
        ax.text(k, ymax + 0.12 * span, st, ha="center", va="bottom",
                fontsize=8)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel(ylabel, fontsize=9)
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)


gene_vals = [(expr.loc[g, s1_ids].astype(float).values,
              expr.loc[g, s2_ids].astype(float).values,
              res.loc[res["gene"] == g, "BH_FDR"].iloc[0]) for g in GENES]
paired_boxes(ax, gene_vals, [DISPLAY[g] for g in GENES],
             "Expression, log2(TPM+1)")
ax.set_xlim(-0.6, len(GENES) - 0.4)
_ymax = max(max(v1.max(), v2.max()) for v1, v2, _ in gene_vals)
_ymax_span = max(max(v1.max(), v2.max()) - min(v1.min(), v2.min())
                 for v1, v2, _ in gene_vals)
ax.set_ylim(top=_ymax + 0.28 * _ymax_span)
for t in ax.get_xticklabels():
    t.set_style("italic")
ax.set_title("Immune/inflammatory PLA2-related genes", fontsize=10)

score_fdr = res.loc[res["gene"] == "SET_SCORE_mean_z", "BH_FDR"].iloc[0]
paired_boxes(ax2, [(s1_sc.values, s2_sc.values, score_fdr)], ["Set\nscore"],
             "Mean z score")
ax2.set_xlim(-0.6, 0.6)
_smax = max(s1_sc.max(), s2_sc.max())
_smin = min(s1_sc.min(), s2_sc.min())
ax2.set_ylim(top=_smax + 0.28 * (_smax - _smin))
for t in ax2.get_xticklabels():
    t.set_style("normal")
ax2.set_title("Score", fontsize=10)

h = [plt.Rectangle((0, 0), 1, 1, fc=C1, alpha=0.55),
     plt.Rectangle((0, 0), 1, 1, fc=C2, alpha=0.55)]
fig.legend(h, [f"S1 (n={len(s1_ids)})", f"S2 (n={len(s2_ids)})"],
           loc="upper center", ncol=2, frameon=False, fontsize=8.5,
           bbox_to_anchor=(0.5, 1.005))
fig.subplots_adjust(top=0.86, bottom=0.28)

os.makedirs(os.path.dirname(FIG_OUT), exist_ok=True)
fig.savefig(FIG_OUT + ".png", dpi=600, bbox_inches="tight")
fig.savefig(FIG_OUT + ".svg", bbox_inches="tight")
print("saved:", TAB_OUT)
print("saved:", FIG_OUT + ".png / .svg")
print(res.round(4).to_string(index=False))