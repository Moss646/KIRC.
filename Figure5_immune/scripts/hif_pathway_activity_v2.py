# hif_pathway_activity_v2.py
# HIF pathway activity between the two lipid-metabolism subtypes
# (S1 = Subtype_1, S2 = Subtype_2) of TCGA-KIRC.
# Layer 1: HALLMARK_HYPOXIA ssGSEA score (Mann-Whitney U).
# Layer 2: expression of 7 HIF genes (Mann-Whitney U, BH across 7).
# Output: results/tables/hif_pathway_activity_stats_v2.csv
#         results/figures/hif_pathway_activity_v2.{png,svg}

import sys, os, io, atexit, string, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)

SCORES_PATH = os.path.join(REPO_ROOT, "data", "processed",
                           "ssgsea_immune_scores_R.csv")
STATS_REF_PATH = os.path.join(REPO_ROOT, "data", "processed",
                              "ssgsea_immune_stats_R.csv")
EXPR_PATH = os.path.join(REPO_ROOT, "data", "tcga_kirc",
                         "KIRC_expr_log2_tpm.csv")
SUB_PATH = os.path.join(REPO_ROOT, "data", "processed",
                        "subtype_assignment_balanced.csv")
OUT_CSV = os.path.join(REPO_ROOT, "results", "tables",
                       "hif_pathway_activity_stats_v2.csv")
OUT_PNG = os.path.join(REPO_ROOT, "results", "figures",
                       "hif_pathway_activity_v2.png")
OUT_SVG = os.path.join(REPO_ROOT, "results", "figures",
                       "hif_pathway_activity_v2.svg")

PATHWAY = "HALLMARK_HYPOXIA"
GENES = ["HIF1A", "EPAS1", "CA9", "VEGFA", "SLC2A1", "EGLN3", "NDUFA4L2"]
S1_COLOR = "#C65A5A"
S2_COLOR = "#4C78A8"
PANELS = [PATHWAY] + GENES

STAR = {0.001: "***", 0.01: "**", 0.05: "*"}


def resolve_gene(name, index):
    """Exact match, then case-insensitive; never silently mismatches."""
    if name in index:
        return name
    lowered = {str(i).lower(): i for i in index}
    hit = lowered.get(name.lower())
    if hit is not None:
        print(f"  resolve_gene: '{name}' -> '{hit}' (case-insensitive)")
    return hit


def mw(s1vals, s2vals):
    return mannwhitneyu(s1vals, s2vals)


def stars(p):
    for th, s in STAR.items():
        if p < th:
            return s
    return "ns"


def main():
    sub = pd.read_csv(SUB_PATH)
    if not {"Patient", "Subtype"}.issubset(sub.columns):
        raise RuntimeError(f"unexpected subtype columns: {list(sub.columns)}")
    lbl = dict(zip(sub["Patient"], sub["Subtype"]))
    s1 = [s for s in sub.loc[sub.Subtype == "Subtype_1", "Patient"]]
    s2 = [s for s in sub.loc[sub.Subtype == "Subtype_2", "Patient"]]
    print(f"samples: S1 n={len(s1)}, S2 n={len(s2)}")

    # ---- Pathway layer ----
    scores = pd.read_csv(SCORES_PATH, index_col=0)
    if PATHWAY not in scores.index:
        raise RuntimeError(f"{PATHWAY} not in {os.path.basename(SCORES_PATH)}")
    hyp = scores.loc[PATHWAY]
    missing_lbl = [s for s in hyp.index if s not in lbl]
    if missing_lbl:
        raise RuntimeError(f"{len(missing_lbl)} score samples lack labels; "
                           f"first: {missing_lbl[:5]}")
    s1p = [s for s in s1 if s in hyp.index]
    s2p = [s for s in s2 if s in hyp.index]
    print(f"  all {len(hyp)} score samples have subtype labels "
          f"(S1 {len(s1p)} / S2 {len(s2p)})")

    stat_ref = pd.read_csv(STATS_REF_PATH)
    pw_col = stat_ref.columns[stat_ref.columns.str.lower()
                              .str.contains("pathway")][0] \
        if stat_ref.columns.str.lower().str.contains("pathway").any() \
        else stat_ref.columns[0]
    row_ref = stat_ref[stat_ref[pw_col] == PATHWAY]
    if len(row_ref) == 0:
        raise RuntimeError(f"{PATHWAY} not in reference stats table")
    row_ref = row_ref.iloc[0]

    rows = []
    u, p_raw = mw(hyp[s1p], hyp[s2p])
    rows.append({
        "layer": "pathway", "item": "HALLMARK_HYPOXIA",
        "S1_mean": round(float(hyp[s1p].mean()), 4),
        "S2_mean": round(float(hyp[s2p].mean()), 4),
        "diff_S1_minus_S2": round(float(hyp[s1p].mean() - hyp[s2p].mean()), 4),
        "p_raw": p_raw,
        "p_bh": float(row_ref["p_bh"]),
        "higher": "S2" if hyp[s1p].mean() < hyp[s2p].mean() else "S1",
    })
    print(f"[pathway] {PATHWAY}: S1={hyp[s1p].mean():.4f} "
          f"S2={hyp[s2p].mean():.4f} p_raw={p_raw:.6g} "
          f"p_bh(ref)={row_ref['p_bh']:.4g}")

    # ---- Gene layer ----
    expr = pd.read_csv(EXPR_PATH, index_col=0)
    unlabelled = [s for s in expr.columns if s not in lbl]
    if unlabelled:
        raise RuntimeError(f"{len(unlabelled)} expression samples lack "
                           f"subtype labels: {unlabelled[:5]}")
    print(f"  all {expr.shape[1]} expression samples have subtype labels")

    genes_ok = []
    for g in GENES:
        hit = resolve_gene(g, expr.index)
        if hit is None:
            print(f"  WARNING: '{g}' not found (exact or case-insensitive); "
                  f"skipped")
        else:
            genes_ok.append(hit)
    if not genes_ok:
        raise RuntimeError("no HIF genes could be resolved; aborting to avoid "
                           "an empty result")

    pvals = [mw(expr.loc[g, s1], expr.loc[g, s2]).pvalue for g in genes_ok]
    p_bh = multipletests(pvals, method="fdr_bh")[1]

    for g, p, q in zip(genes_ok, pvals, p_bh):
        m1, m2 = expr.loc[g, s1].mean(), expr.loc[g, s2].mean()
        rows.append({
            "layer": "gene", "item": g,
            "S1_mean": round(float(m1), 4),
            "S2_mean": round(float(m2), 4),
            "diff_S1_minus_S2": round(float(m1 - m2), 4),
            "p_raw": p,
            "p_bh": q,
            "higher": "S2" if m1 < m2 else "S1",
        })
        print(f"[gene] {g}: S1={m1:.4f} S2={m2:.4f} "
              f"p_raw={p:.3g} p_bh={q:.3g}")

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote: {OUT_CSV}")

    # ---- Figure: adaptive boxplot grid ----
    n = len(PANELS)
    ncols = 4
    nrows = int(math.ceil(n / ncols))
    letters = string.ascii_uppercase[:n]
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols,
                                                    3.2 * nrows))
    axes = np.array(axes).ravel()
    for ax in axes[n:]:
        ax.set_visible(False)
    for i, (ax, item) in enumerate(zip(axes, PANELS)):
        if item == PATHWAY:
            vals1, vals2 = hyp[s1p], hyp[s2p]
            ylab = "ssGSEA score"
            title = "HALLMARK\nHYPOXIA"
        else:
            vals1, vals2 = expr.loc[item, s1], expr.loc[item, s2]
            ylab = "log2(TPM+1)"
            title = item
        bp = ax.boxplot([vals1, vals2], positions=[0, 1], widths=0.55,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.4))
        for patch, c in zip(bp["boxes"], [S1_COLOR, S2_COLOR]):
            patch.set_facecolor(c)
            patch.set_alpha(0.65)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["S1", "S2"], fontsize=10)
        ax.set_title(title, fontsize=11, pad=8)
        ax.text(-0.14, 1.02, letters[i], transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="bottom", ha="left")
        ax.set_ylabel(ylab, fontsize=10)
        u, p = mw(vals1, vals2)
        ymax = max(np.percentile(vals1, 97), np.percentile(vals2, 97))
        ax.plot([0, 0, 1, 1], [ymax, ymax * 1.06, ymax * 1.06, ymax],
                color="black", lw=1.0)
        ax.text(0.5, ymax * 1.10, stars(p), ha="center",
                va="bottom", fontsize=12)
        ax.tick_params(labelsize=9)

    fig.suptitle("HIF pathway activity: S1 vs S2", fontsize=13, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=600)
    fig.savefig(OUT_SVG)
    plt.close(fig)
    print(f"wrote: {OUT_PNG}")
    print(f"wrote: {OUT_SVG}")


if __name__ == "__main__":
    main()