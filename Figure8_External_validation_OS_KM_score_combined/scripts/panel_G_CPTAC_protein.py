"""CPTAC protein validation of the 14-gene lipid metabolism panel.

Outputs:
    results/figures/Figure_G_CPTAC_protein_v1.png
    results/figures/Figure_G_CPTAC_protein_v1.svg
"""

import io
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"

from config import *

FIG = FIGURES_DIR
C1 = "#C65A5A"
C2 = "#4C78A8"

PANEL14 = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36",
    "SLC27A2", "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]


def stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def main():
    rnaseq_cptac = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
    prot_tumor = pd.read_csv(CPTAC_PROTEIN, sep="\t", index_col=0)

    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    avail_clf_genes = [g for g in clf["genes"] if g in rnaseq_cptac.index]
    avail_clf_idx = [clf["genes"].index(g) for g in avail_clf_genes]

    expr_clf = rnaseq_cptac.loc[avail_clf_genes].values.T
    expr_clf_z = (expr_clf - expr_clf.mean(axis=0)) / (
        expr_clf.std(axis=0, ddof=1) + 1e-8
    )

    c1 = np.array(clf["centroid_s1"])[avail_clf_idx]
    c2 = np.array(clf["centroid_s2"])[avail_clf_idx]
    d1 = cdist(expr_clf_z, c1.reshape(1, -1))[:, 0]
    d2 = cdist(expr_clf_z, c2.reshape(1, -1))[:, 0]
    cptac_labels = np.where(d1 < d2, 1, 0)

    cptac_samps = list(rnaseq_cptac.columns)
    cptac_s1 = np.array(cptac_samps)[cptac_labels == 1]
    cptac_s2 = np.array(cptac_samps)[cptac_labels == 0]

    print(
        f"CPTAC NCC: S1={len(cptac_s1)}, S2={len(cptac_s2)}, "
        f"genes={len(avail_clf_genes)}/68"
    )

    f_avail = [g for g in PANEL14 if g in prot_tumor.index]
    print(f"Protein genes available: {len(f_avail)}/14")

    rows = []
    for g in f_avail:
        s1 = prot_tumor.loc[
            g, [c for c in cptac_s1 if c in prot_tumor.columns]
        ].values.astype(float)
        s2 = prot_tumor.loc[
            g, [c for c in cptac_s2 if c in prot_tumor.columns]
        ].values.astype(float)

        s1 = s1[~np.isnan(s1)]
        s2 = s2[~np.isnan(s2)]

        if len(s1) < 3 or len(s2) < 3:
            print(f"  {g:10s}: skipped (S1={len(s1)}, S2={len(s2)})")
            continue

        fc = s1.mean() - s2.mean()
        se = np.sqrt(np.var(s1) / len(s1) + np.var(s2) / len(s2))
        _, p = mannwhitneyu(s1, s2)

        rows.append({
            "gene": g,
            "log2FC": fc,
            "ci_low": fc - 1.96 * se,
            "ci_high": fc + 1.96 * se,
            "p": p,
            "S1_mean": s1.mean(),
            "S2_mean": s2.mean(),
            "n_S1": len(s1),
            "n_S2": len(s2),
        })

        direction = "S1>" if fc > 0 else "S2>"
        print(
            f"  {g:10s}: S1={s1.mean():+.4f} S2={s2.mean():+.4f} "
            f"FC={fc:+.4f} P={p:.2e} {stars(p)} [{direction}] "
            f"n=({len(s1)},{len(s2)})"
        )

    f_df = pd.DataFrame(rows).sort_values("log2FC")

    if len(f_df) > 0:
        _, f_df["p_adj"], _, _ = multipletests(f_df["p"].values, method="fdr_bh")
        f_df["sig_adj"] = f_df["p_adj"].apply(stars)
        n_bh = int((f_df["p_adj"] < 0.05).sum())
        print(f"\nBH FDR<0.05: {n_bh}/{len(f_df)}")
        for _, r in f_df.iterrows():
            print(
                f"  {r['gene']:10s}: FC={r['log2FC']:+.4f} "
                f"raw_P={r['p']:.2e} BH_P={r['p_adj']:.2e} {r['sig_adj']}"
            )

    texpr = pd.read_csv(TCGA_EXPR, index_col=0)
    tsub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)
    t1 = [c for c in tsub[tsub.Subtype == "Subtype_1"]["Patient"] if c in texpr.columns]
    t2 = [c for c in tsub[tsub.Subtype == "Subtype_2"]["Patient"] if c in texpr.columns]

    tcga_fc = {}
    for g in PANEL14:
        if g in texpr.index:
            tcga_fc[g] = (
                texpr.loc[g, t1].astype(float).mean()
                - texpr.loc[g, t2].astype(float).mean()
            )

    valid_genes = [g for g in f_df["gene"] if g in tcga_fc]
    agreement = sum(
        (f_df.loc[f_df.gene == g, "log2FC"].values[0] > 0) == (tcga_fc[g] > 0)
        for g in valid_genes
    )
    print(f"\nDirection agreement with TCGA: {agreement}/{len(valid_genes)}")

    fig = plt.figure(figsize=(10, 8))
    fig.patch.set_facecolor("white")
    ax = plt.subplot(111)
    ax.set_facecolor("white")

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    f_plot = f_df.sort_values("log2FC")
    genes_g = list(f_plot["gene"])[::-1]
    fcs_g = list(f_plot["log2FC"])[::-1]
    los_g = list(f_plot["ci_low"])[::-1]
    his_g = list(f_plot["ci_high"])[::-1]
    y_g = np.arange(len(genes_g))

    for i, g in enumerate(genes_g):
        col = C1 if fcs_g[i] > 0 else C2
        ax.barh(i, fcs_g[i], color=col, alpha=0.85, height=0.55, zorder=3)

    ax.errorbar(
        fcs_g,
        y_g,
        xerr=[
            [fc - cl for fc, cl in zip(fcs_g, los_g)],
            [ch - fc for fc, ch in zip(fcs_g, his_g)],
        ],
        fmt="none",
        ecolor="#333",
        capsize=2,
        lw=0.8,
        alpha=0.6,
        zorder=4,
    )

    ax.set_yticks(y_g)
    ax.set_yticklabels(genes_g, fontsize=14, fontweight="bold", fontstyle="italic")
    ax.axvline(x=0, color="#333", lw=1.2)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xlabel(r"CPTAC Protein log$_2$ FC (S1 minus S2)", fontsize=20)
    ax.set_title(
        "CPTAC Proteomic Validation of Lipid Metabolism Genes",
        fontsize=16,
        fontweight="bold",
        pad=10,
    )
    ax.legend(
        handles=[Patch(facecolor=C1, alpha=0.85), Patch(facecolor=C2, alpha=0.85)],
        labels=["S1", "S2"],
        fontsize=14,
        loc="upper right",
        framealpha=0.9,
    )
    ax.tick_params(labelsize=18)

    out_png = os.path.join(FIG, "Figure_G_CPTAC_protein_v1.png")
    out_svg = os.path.join(FIG, "Figure_G_CPTAC_protein_v1.svg")
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"\nSaved: {os.path.basename(out_png)}")
    print(f"       {os.path.basename(out_svg)}")

    try:
        n_valid = len(valid_genes)
        n_bh = int((f_df["p_adj"] < 0.05).sum())
        n_tot = len(f_df)
        prov = [
            {
                "section": "Biological validation",
                "item": "14-gene direction concordance",
                "cohort": "CPTAC",
                "value": f"{agreement} of {n_valid}",
            },
            {
                "section": "Biological validation",
                "item": "14-gene BH-sig genes (MWU FDR<0.05)",
                "cohort": "CPTAC",
                "value": f"{n_bh} of {n_tot}",
            },
        ]
        emit_provenance("panel_G_CPTAC_protein.py", prov)
    except Exception as e:
        print(f"[PROVENANCE] panel_G_CPTAC_protein emission skipped: {e}")


if __name__ == "__main__":
    main()