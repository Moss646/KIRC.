"""9-gene non-NCC FAO panel: four-cohort S1/S2 analysis.

Outputs:
    results/tables/nonncc9_genes_comparison_v1.csv
    results/tables/nonncc9_cross_cohort_summary_v1.csv
    results/tables/nonncc9_fao_score_stats_v1.csv
    results/tables/nonncc9_fao_score_per_sample_v1.csv
    results/figures/nonncc9_fao_score_4cohorts_v1.png / .svg
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, ttest_ind
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from config import (
    TABLES_DIR,
    FIGURES_DIR,
    CPTAC_PROTEIN,
    CPTAC_RNA,
    SUBTYPE_CLASSIFIER_JSON,
    SUBTYPE_ASSIGNMENT_CSV,
    TCGA_EXPR,
    TCGA_DEG_CSV,
)
from key14_cross_cohort_v2 import (
    assign_subtypes,
    load_mapping,
    prepare_emtab,
    prepare_icgc,
)

NON_NCC9 = [
    "CPT1A", "ACOX1", "CPT2", "CD36", "SLC27A2",
    "FASN", "SCD", "PLIN2", "PPARG",
]

RELAY = Path(os.environ.get(
    "NONNCC9_RELAY",
    Path(tempfile.gettempdir()) / "nonncc9",
))

R_SCRIPT = SCRIPT_DIR / "nonncc9_limma_external_v1.R"

RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript")
if not RSCRIPT:
    raise SystemExit(
        "Rscript not found. Install R, add it to PATH, "
        "or set the RSCRIPT environment variable."
    )

C1 = "#C65A5A"
C2 = "#4C78A8"


def tcga_reference():
    limma = pd.read_csv(TCGA_DEG_CSV)
    li = limma.set_index("gene_symbol")
    return {
        g: (float(li.loc[g, "logFC"]), float(li.loc[g, "FDR"]))
        for g in NON_NCC9
        if g in li.index
    }


def zscore_matrix(mat):
    mu = mat.mean(axis=0)
    sd = mat.std(axis=0, ddof=1)
    return mat.subtract(mu, axis=1).div(sd.replace(0, np.nan), axis=1)


def cohens_d(v1, v2):
    n1, n2 = len(v1), len(v2)
    s = np.sqrt(
        ((n1 - 1) * v1.std(ddof=1) ** 2 + (n2 - 1) * v2.std(ddof=1) ** 2)
        / (n1 + n2 - 2)
    )
    return (v1.mean() - v2.mean()) / s if s > 0 else np.nan


def main():
    RELAY.mkdir(parents=True, exist_ok=True)
    Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)
    Path(TABLES_DIR).mkdir(parents=True, exist_ok=True)

    for p in [TCGA_DEG_CSV, TCGA_EXPR, SUBTYPE_ASSIGNMENT_CSV,
              SUBTYPE_CLASSIFIER_JSON, CPTAC_PROTEIN, CPTAC_RNA, R_SCRIPT]:
        if not Path(p).exists():
            raise SystemExit(f"missing input: {p}")

    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    ref_genes = clf["genes"]
    tref = tcga_reference()
    print(f"TCGA reference: {len(tref)}/9 non-NCC genes in limma table")

    # TCGA score
    tpm = pd.read_csv(TCGA_EXPR, index_col=0)
    labels = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)
    labels["group"] = labels["Subtype"].map(
        {"Subtype_1": "S1", "Subtype_2": "S2"}
    )
    lab = labels.set_index("Patient")["group"]

    common = [c for c in tpm.columns if c in lab.index]
    tsub = tpm.loc[[g for g in NON_NCC9 if g in tpm.index], common]
    tz = zscore_matrix(tsub.T)
    t_score = tz.mean(axis=1)
    t_grp = lab.loc[t_score.index]

    # E-MTAB
    g2e = load_mapping()
    emat = prepare_emtab()
    e_s1i, e_s2i = assign_subtypes(emat, ref_genes, clf)
    print(f"E-MTAB: {emat.shape[0]} samples, "
          f"S1-like={len(e_s1i)} S2-like={len(e_s2i)}")

    e_genes = [g for g in NON_NCC9 if g in emat.columns]
    e_expr = emat[e_genes].iloc[e_s1i + e_s2i]
    e_grpv = (["S1like"] * len(e_s1i)) + (["S2like"] * len(e_s2i))
    e_score = zscore_matrix(emat[e_genes]).iloc[e_s1i + e_s2i].mean(axis=1)
    e_score.index = e_expr.index
    e_grp = pd.Series(e_grpv, index=e_expr.index)

    # ICGC
    all_sym = list(set(NON_NCC9) | set(ref_genes))
    iem = prepare_icgc(g2e, all_sym)
    i_s1i, i_s2i = assign_subtypes(iem, ref_genes, clf)
    print(f"ICGC:   {iem.shape[0]} samples, "
          f"S1-like={len(i_s1i)} S2-like={len(i_s2i)}")

    i_genes = [g for g in NON_NCC9 if g in iem.columns]
    i_expr = iem[i_genes].iloc[i_s1i + i_s2i]
    i_grpv = (["S1like"] * len(i_s1i)) + (["S2like"] * len(i_s2i))
    i_score = zscore_matrix(iem[i_genes]).iloc[i_s1i + i_s2i].mean(axis=1)
    i_score.index = i_expr.index
    i_grp = pd.Series(i_grpv, index=i_expr.index)

    # CPTAC protein
    prot = pd.read_csv(CPTAC_PROTEIN, sep="\t", index_col=0)
    rnaseq_cp = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
    rnaz = rnaseq_cp.subtract(rnaseq_cp.mean(axis=1), axis=0).div(
        rnaseq_cp.std(axis=1, ddof=1), axis=0
    ).fillna(0)

    cp_s1i, cp_s2i = assign_subtypes(rnaz.T, ref_genes, clf)
    cp_s1 = rnaseq_cp.columns[cp_s1i]
    cp_s2 = rnaseq_cp.columns[cp_s2i]
    print(f"CPTAC:  protein S1-like={len(cp_s1)} S2-like={len(cp_s2)}")

    cp_cols = list(cp_s1) + list(cp_s2)
    cp_genes = [g for g in NON_NCC9 if g in prot.index]
    cp_expr = prot.loc[cp_genes, [c for c in cp_cols if c in prot.columns]].T
    cp_grpv = (["S1like"] * len(cp_s1)) + (["S2like"] * len(cp_s2))
    cp_grp = pd.Series(cp_grpv, index=cp_expr.index)
    cp_score = zscore_matrix(cp_expr).mean(axis=1)

    # Export for R limma
    e_expr.T.to_csv(RELAY / "emtab_expr.csv")
    e_grp.rename("group").to_frame().assign(sample=e_grp.index).to_csv(
        RELAY / "emtab_groups.csv", index=False
    )
    i_expr.T.to_csv(RELAY / "icgc_expr.csv")
    i_grp.rename("group").to_frame().assign(sample=i_grp.index).to_csv(
        RELAY / "icgc_groups.csv", index=False
    )
    cp_expr.T.to_csv(RELAY / "cptac_expr.csv")
    cp_grp.rename("group").to_frame().assign(sample=cp_grp.index).to_csv(
        RELAY / "cptac_groups.csv", index=False
    )

    r_local = RELAY / "nonncc9_limma_external_v1.R"
    shutil.copy(R_SCRIPT, r_local)

    print("Running limma (R) ...")
    r = subprocess.run(
        [RSCRIPT, str(r_local), str(RELAY)],
        capture_output=True,
        text=True,
    )
    print(r.stdout[-1500:])
    if r.returncode != 0:
        print(r.stderr[-3000:])
        raise SystemExit("R limma failed")

    rows = []
    for g in NON_NCC9:
        lf, fdr = tref[g]
        rows.append({
            "gene": g,
            "cohort": "TCGA",
            "diff": lf,
            "p": np.nan,
            "p_adj": fdr,
            "effect_measure": "limma_voom_logFC (reference)",
            "dir_match": "",
        })

    lim = pd.read_csv(RELAY / "nonncc9_limma_external_v1.csv")
    for c in ["E-MTAB", "ICGC", "CPTAC"]:
        sub = lim[lim.cohort == c]
        for _, rrow in sub.iterrows():
            g = rrow["gene"]
            rows.append({
                "gene": g,
                "cohort": c,
                "diff": rrow["logFC"],
                "p": rrow["P.Value"],
                "p_adj": rrow["adj.P.Val"],
                "effect_measure": "limma_logFC",
                "dir_match": int((rrow["logFC"] > 0) == (tref[g][0] > 0)),
            })

    comp = pd.DataFrame(rows)
    comp.to_csv(Path(TABLES_DIR) / "nonncc9_genes_comparison_v1.csv", index=False)

    tcga_dir = comp[comp.cohort == "TCGA"].set_index("gene")["diff"]
    sumrows = []

    for c in ["TCGA", "E-MTAB", "ICGC", "CPTAC"]:
        sub = comp[comp.cohort == c].dropna(subset=["diff"])
        n_det = len(sub)

        if c == "TCGA":
            nsig = int((sub.p_adj < 0.05).sum())
            sumrows.append({
                "cohort": c, "n_detected": n_det, "dir_match": "",
                "n_sig": nsig, "sig_dir_match": "",
                "pearson_n": "", "pearson_r": "", "pearson_P": "",
            })
            continue

        match = int((sub.dir_match == 1).sum())
        sig = sub[sub.p_adj < 0.05]
        sig_match = int((sig.dir_match == 1).sum())

        ext = sub.set_index("gene")["diff"]
        shared = ext.index.intersection(tcga_dir.index)
        rr = pearsonr(ext.loc[shared], tcga_dir.loc[shared])

        sumrows.append({
            "cohort": c,
            "n_detected": n_det,
            "dir_match": f"{match}/{n_det}",
            "n_sig": len(sig),
            "sig_dir_match": f"{sig_match}/{len(sig)}",
            "pearson_n": len(shared),
            "pearson_r": round(rr[0], 3),
            "pearson_P": f"{rr[1]:.2e}",
        })

    summ = pd.DataFrame(sumrows)
    summ.to_csv(Path(TABLES_DIR) / "nonncc9_cross_cohort_summary_v1.csv", index=False)
    print("\nPer-cohort summary")
    print(summ.to_string(index=False))

    def score_stats(name, score, grp, n_genes):
        v1 = score[grp.isin(["S1", "S1like"])].dropna()
        v2 = score[grp.isin(["S2", "S2like"])].dropna()
        mw = mannwhitneyu(v1, v2, alternative="two-sided")[1]
        tt = ttest_ind(v1, v2, equal_var=False)[1]
        return {
            "cohort": name,
            "n_genes_used": n_genes,
            "n_s1": len(v1),
            "n_s2": len(v2),
            "mean_S1": v1.mean(),
            "sd_S1": v1.std(ddof=1),
            "mean_S2": v2.mean(),
            "sd_S2": v2.std(ddof=1),
            "diff_S1_minus_S2": v1.mean() - v2.mean(),
            "cohens_d": cohens_d(v1, v2),
            "p_mwu": mw,
            "p_welch": tt,
        }

    stats = pd.DataFrame([
        score_stats("TCGA", t_score, t_grp, 9),
        score_stats("E-MTAB", e_score, e_grp, len(e_genes)),
        score_stats("ICGC", i_score, i_grp, len(i_genes)),
        score_stats("CPTAC", cp_score, cp_grp, len(cp_genes)),
    ])
    stats["p_mwu_BH_4cohorts"] = multipletests(
        stats.p_mwu, method="fdr_bh"
    )[1]
    stats.to_csv(Path(TABLES_DIR) / "nonncc9_fao_score_stats_v1.csv", index=False)
    print("\n9-gene non-NCC FAO score")
    print(stats.round(4).to_string(index=False))

    per_sample = []
    for name, sc, gr in [
        ("TCGA", t_score, t_grp),
        ("E-MTAB", e_score, e_grp),
        ("ICGC", i_score, i_grp),
        ("CPTAC", cp_score, cp_grp),
    ]:
        per_sample.append(
            pd.DataFrame({
                "cohort": name,
                "sample": sc.index,
                "group": gr.values,
                "score": sc.values,
            })
        )
    pd.concat(per_sample).to_csv(
        Path(TABLES_DIR) / "nonncc9_fao_score_per_sample_v1.csv", index=False
    )

    plt.rcParams.update({
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "axes.linewidth": 0.8,
    })

    panels = [
        ("TCGA", "TCGA", t_score, t_grp,
         f"S1 n={int((t_grp == 'S1').sum())}\nS2 n={int((t_grp == 'S2').sum())}"),
        ("E-MTAB", "E-MTAB-1980", e_score, e_grp,
         f"S1-like n={int((e_grp == 'S1like').sum())}\nS2-like n={int((e_grp == 'S2like').sum())}"),
        ("ICGC", "ICGC RECA-EU", i_score, i_grp,
         f"S1-like n={int((i_grp == 'S1like').sum())}\nS2-like n={int((i_grp == 'S2like').sum())}"),
        ("CPTAC", "CPTAC protein", cp_score, cp_grp,
         f"S1-like n={int((cp_grp == 'S1like').sum())}\nS2-like n={int((cp_grp == 'S2like').sum())}"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(7.087, 2.9))

    for ax, (key, name, sc, gr, ntxt) in zip(axes, panels):
        v1 = sc[gr.isin(["S1", "S1like"])].dropna()
        v2 = sc[gr.isin(["S2", "S2like"])].dropna()

        for vals, pos, col in [(v1, 0, C1), (v2, 1, C2)]:
            ax.boxplot(
                vals,
                positions=[pos],
                widths=0.55,
                patch_artist=True,
                medianprops=dict(color="black", lw=1.0),
                boxprops=dict(facecolor=col, alpha=0.55, lw=0.8),
                whiskerprops=dict(lw=0.8),
                capprops=dict(lw=0.8),
                flierprops=dict(marker="", lw=0),
            )
            rng = np.random.default_rng(42 + pos)
            ax.scatter(
                rng.uniform(pos - 0.18, pos + 0.18, len(vals)),
                vals,
                s=3,
                color=col,
                alpha=0.45,
                lw=0,
            )

        p = stats.loc[stats.cohort == key, "p_mwu_BH_4cohorts"].iloc[0]
        top = max(v1.max(), v2.max())
        ax.set_ylim(top=top + 0.55)
        ax.text(0.5, top + 0.34, f"P_BH = {p:.1e}", ha="center", fontsize=8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(
            ["S1", "S2"] if name == "TCGA" else ["S1-like", "S2-like"]
        )
        ax.set_title(f"{name}\n{ntxt}", fontsize=10)
        ax.set_xlim(-0.6, 1.6)
        ax.tick_params(labelsize=9, width=0.8)

        for sp in ax.spines.values():
            sp.set_linewidth(0.8)

    axes[0].set_ylabel(
        "9-gene non-NCC FAO score\n(mean z; higher = S2-like)", fontsize=10
    )
    fig.tight_layout()
    fig.savefig(Path(FIGURES_DIR) / "nonncc9_fao_score_4cohorts_v1.png", dpi=600)
    fig.savefig(Path(FIGURES_DIR) / "nonncc9_fao_score_4cohorts_v1.svg")
    plt.close(fig)
    print("\nDone.")


if __name__ == "__main__":
    main()