"""14-gene lipid-metabolism panel: cross-cohort validation.

Uniform mean-difference reference across cohorts. TCGA reference is not on
limma, matching key14_biological.py and lipid16_cross_cohort_v2.py.

Effect estimator (same for every cohort):
    S1 mean minus S2 mean on the cohort's native log scale.

Cohorts:
    TCGA-KIRC    log2TPM, frozen subtype assignment, MWU + BH
    E-MTAB-1980  microarray log intensities, nearest-centroid, MWU + BH
    ICGC RECA-EU log2(normalized_read_count + 1), same assignment, MWU + BH
    CPTAC        protein log abundance, effect only (no counts, no MWU here)

Pearson r: external effect vs TCGA mean difference. Fisher-z 95% CI.
Sign convention: >0 means higher in S1.

Outputs:
    results/tables/key14_genes_comparison_v2.csv
    results/tables/key14_cross_cohort_summary_v2.csv
"""

import gzip
import json
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, pearsonr
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from config import (
    GEO_EMTAB_EXPR,
    GEO_NM2GENE,
    GENE_MAPPING_JSON,
    ICGC_EXP_SEQ,
    ICGC_SPECIMEN,
    CPTAC_PROTEIN,
    CPTAC_RNA,
    SUBTYPE_CLASSIFIER_JSON,
    TCGA_EXPR,
    SUBTYPE_ASSIGNMENT_CSV,
    TABLES_DIR,
    emit_provenance,
)

KEY14 = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36",
    "SLC27A2", "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]


def load_mapping():
    with open(GENE_MAPPING_JSON) as f:
        return json.load(f)


def prepare_emtab():
    expr = pd.read_csv(GEO_EMTAB_EXPR, sep="\t", index_col=0)
    sys_col = expr["SystematicName"].values

    with open(GEO_NM2GENE) as f:
        nm2g = json.load(f)

    gene_data = defaultdict(list)
    for i in range(len(sys_col)):
        nm = sys_col[i]
        if nm.startswith("NM_") and nm in nm2g:
            gene_data[nm2g[nm]].append(expr.iloc[i, 2:].astype(float).values)

    return pd.DataFrame(
        {g: np.mean(arrs, axis=0) for g, arrs in gene_data.items()},
        index=[c.split("_")[0] for c in expr.columns[2:]],
    )


def prepare_icgc(g2e, symbols):
    e2s = {v: k for k, v in g2e.items()}

    tumor_dons = set()
    with gzip.open(ICGC_SPECIMEN, "rt") as f:
        hdr = f.readline().strip().split("\t")
        sc = hdr.index("specimen_type")
        did_idx = hdr.index("icgc_donor_id")
        for line in f:
            c = line.strip().split("\t")
            if "Primary tumour" in c[sc]:
                tumor_dons.add(c[did_idx])

    iexp = defaultdict(lambda: defaultdict(float))
    with gzip.open(ICGC_EXP_SEQ, "rt") as f:
        hdr = f.readline().strip().split("\t")
        gi = hdr.index("gene_id")
        si = hdr.index("icgc_donor_id")
        fi = hdr.index("normalized_read_count")
        for line in f:
            c = line.strip().split("\t")
            gv = c[gi].replace("gene:", "")
            if gv in e2s and c[si] in tumor_dons:
                iexp[c[si]][e2s[gv]] = float(c[fi])

    iem = pd.DataFrame(
        {d: {g: iexp[d].get(g, 0) for g in symbols if g in g2e} for d in iexp}
    ).T.astype(float)

    return np.log2(iem + 1)


def assign_subtypes(mat, ref_genes, clf):
    mz = mat.subtract(mat.mean(axis=0), axis=1).div(
        mat.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    c1 = np.array(clf["centroid_s1"])
    c2 = np.array(clf["centroid_s2"])
    com = [g for g in ref_genes if g in mz.columns]
    tci = [ref_genes.index(g) for g in com]

    d1 = cdist(mz[com].values, c1[tci].reshape(1, -1))[:, 0]
    d2 = cdist(mz[com].values, c2[tci].reshape(1, -1))[:, 0]

    return list(np.where(d1 < d2)[0]), list(np.where(d1 >= d2)[0])


def fc_map(mat, s1i, s2i, genes):
    out = {}
    for g in genes:
        if g not in mat.columns:
            continue
        v1 = mat.iloc[s1i][g].astype(float).dropna()
        v2 = mat.iloc[s2i][g].astype(float).dropna()
        if len(v1) >= 3 and len(v2) >= 3:
            out[g] = v1.mean() - v2.mean()
    return out


def tcga_mean_diff():
    expr = pd.read_csv(TCGA_EXPR, index_col=0)
    sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)

    s1 = [c for c in sub[sub.Subtype == "Subtype_1"]["Patient"] if c in expr.columns]
    s2 = [c for c in sub[sub.Subtype == "Subtype_2"]["Patient"] if c in expr.columns]

    out = {}
    for g in KEY14:
        if g in expr.index:
            v1 = expr.loc[g, s1].astype(float)
            v2 = expr.loc[g, s2].astype(float)
            out[g] = (float(v1.mean() - v2.mean()), float(mannwhitneyu(v1, v2)[1]))

    return out


def main():
    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    ref_genes = clf["genes"]

    tmd = tcga_mean_diff()
    print(f"TCGA: {len(tmd)}/14 genes in log2TPM matrix")

    g2e = load_mapping()

    emat = prepare_emtab()
    e_s1i, e_s2i = assign_subtypes(emat, ref_genes, clf)
    print(f"E-MTAB: {emat.shape[0]} samples, S1={len(e_s1i)} S2={len(e_s2i)}")
    em_fc = fc_map(emat, e_s1i, e_s2i, KEY14)

    all_sym = list(set(KEY14) | set(ref_genes))
    iem = prepare_icgc(g2e, all_sym)
    i_s1i, i_s2i = assign_subtypes(iem, ref_genes, clf)
    print(f"ICGC:   {iem.shape[0]} samples, S1={len(i_s1i)} S2={len(i_s2i)}")
    ic_fc = fc_map(iem, i_s1i, i_s2i, KEY14)

    prot = pd.read_csv(CPTAC_PROTEIN, sep="\t", index_col=0)
    rnaseq_cp = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
    rnaz = rnaseq_cp.subtract(rnaseq_cp.mean(axis=1), axis=0).div(
        rnaseq_cp.std(axis=1, ddof=1), axis=0
    ).fillna(0)

    cp_s1i, cp_s2i = assign_subtypes(rnaz.T, ref_genes, clf)
    cp_s1 = rnaseq_cp.columns[cp_s1i]
    cp_s2 = rnaseq_cp.columns[cp_s2i]
    print(f"CPTAC:  protein S1={len(cp_s1)} S2={len(cp_s2)}")

    cp_fc = {}
    for g in KEY14:
        if g not in prot.index:
            continue
        v1 = prot.loc[g, [c for c in cp_s1 if c in prot.columns]].astype(float).dropna()
        v2 = prot.loc[g, [c for c in cp_s2 if c in prot.columns]].astype(float).dropna()
        if len(v1) >= 3 and len(v2) >= 3:
            cp_fc[g] = v1.mean() - v2.mean()

    rows = []
    for g in KEY14:
        diff, p = tmd.get(g, (np.nan, np.nan))
        rows.append({
            "gene": g,
            "cohort": "TCGA",
            "diff": diff,
            "p": p,
            "effect_measure": "mean_diff_log2TPM",
        })

    for cohort, fc in [("E-MTAB", em_fc), ("ICGC", ic_fc), ("CPTAC", cp_fc)]:
        for g in KEY14:
            rows.append({
                "gene": g,
                "cohort": cohort,
                "diff": fc.get(g, np.nan),
                "p_adj": np.nan,
                "effect_measure": "mean_diff_log_scale",
            })

    rs = pd.DataFrame(rows)

    # MWU + BH within E-MTAB and ICGC. CPTAC protein p is not computed here.
    for cohort, mat, s1i, s2i in [
        ("E-MTAB", emat, e_s1i, e_s2i),
        ("ICGC", iem, i_s1i, i_s2i),
    ]:
        for g in KEY14:
            try:
                v1 = mat.iloc[s1i][g].astype(float)
                v2 = mat.iloc[s2i][g].astype(float)
                p = mannwhitneyu(v1, v2)[1]
            except KeyError:
                p = np.nan

            m = (rs["cohort"] == cohort) & (rs["gene"] == g)
            rs.loc[m, "p"] = p

        mask = rs["cohort"] == cohort
        pv = rs.loc[mask, "p"].dropna()
        if len(pv) > 0:
            _, padj, _, _ = multipletests(pv, method="fdr_bh")
            rs.loc[pv.index, "p_adj"] = padj

    mask = rs["cohort"] == "TCGA"
    pv = rs.loc[mask, "p"].dropna()
    if len(pv) > 0:
        _, padj, _, _ = multipletests(pv, method="fdr_bh")
        rs.loc[pv.index, "p_adj"] = padj

    print("\nDirection concordance (sign of S1-S2 vs TCGA mean diff)")
    summary = []
    tcga_dir = rs[rs.cohort == "TCGA"].set_index("gene")["diff"]

    for c in ["TCGA", "E-MTAB", "ICGC", "CPTAC"]:
        sub = rs[rs.cohort == c].dropna(subset=["diff"])
        avail = len(sub)

        if c == "TCGA":
            sig = sub[sub.p_adj < 0.05]
            print(f"{c:6s}: {avail}/14 genes, MWU BH<0.05: {len(sig)}")
            summary.append({
                "cohort": c,
                "n_detected": avail,
                "dir_match": np.nan,
                "n_bh_sig": len(sig),
                "bh_dir_match": np.nan,
            })
            continue

        match = int(
            ((sub["diff"] > 0).values == (tcga_dir.reindex(sub["gene"]) > 0).values).sum()
        )

        if c in ("E-MTAB", "ICGC"):
            bh = sub[sub.p_adj < 0.05]
        else:
            bh = sub

        if len(bh):
            bh_match = int(
                ((bh["diff"] > 0).values == (tcga_dir.reindex(bh["gene"]) > 0).values).sum()
            )
        else:
            bh_match = 0

        print(
            f"{c:6s}: {avail}/14 genes detected, direction match {match}/{avail}, "
            f"BH-sig {len(bh)}/14 (direction match {bh_match}/{len(bh)})"
        )
        summary.append({
            "cohort": c,
            "n_detected": avail,
            "dir_match": match,
            "n_bh_sig": len(bh),
            "bh_dir_match": bh_match,
        })

    print("\nEffect-size Pearson r vs TCGA mean difference")
    for c in ["E-MTAB", "ICGC", "CPTAC"]:
        sub = rs[rs.cohort == c].dropna(subset=["diff"])
        common = sorted(set(sub.gene) & set(tcga_dir.index))
        x = [tcga_dir[g] for g in common]
        y = [sub.set_index("gene").loc[g, "diff"] for g in common]

        r, p = pearsonr(x, y)
        z = np.arctanh(r)
        se = 1 / np.sqrt(len(common) - 3)
        lo = np.tanh(z - 1.96 * se)
        hi = np.tanh(z + 1.96 * se)

        print(f"{c:6s}: n={len(common)} r={r:.3f} ({lo:.3f}-{hi:.3f}) P={p:.1e}")
        summary.append({
            "cohort": c,
            "pearson_n": len(common),
            "pearson_r": r,
            "pearson_CI": f"{lo:.3f}-{hi:.3f}",
            "pearson_P": p,
        })

    out = os.path.join(TABLES_DIR, "key14_genes_comparison_v2.csv")
    rs.to_csv(out, index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(TABLES_DIR, "key14_cross_cohort_summary_v2.csv"), index=False)
    print(f"\nSaved: {out}")
    print(f"Saved: {os.path.join(TABLES_DIR, 'key14_cross_cohort_summary_v2.csv')}")

    prov = []
    cname = {"E-MTAB": "E-MTAB-1980", "ICGC": "ICGC RECA-EU", "CPTAC": "CPTAC"}

    for s in summary:
        if s["cohort"] not in cname:
            continue
        c = cname[s["cohort"]]

        if s.get("dir_match") is not None:
            prov.append({
                "section": "Biological validation",
                "item": "14-gene direction concordance",
                "cohort": c,
                "value": f"{s['dir_match']} of {s['n_detected']}",
            })

        if s.get("pearson_r") is not None:
            prov.append({
                "section": "Biological validation",
                "item": "14-gene Pearson r vs TCGA",
                "cohort": c,
                "value": f"+{s['pearson_r']:.3f}",
            })
            prov.append({
                "section": "Biological validation",
                "item": "14-gene 95% CI",
                "cohort": c,
                "value": s["pearson_CI"],
            })
            prov.append({
                "section": "Biological validation",
                "item": "14-gene Pearson P",
                "cohort": c,
                "value": f"{s['pearson_P']:.1e}",
            })

    emit_provenance("key14_cross_cohort_v2.py", prov)

    print("\nPer-gene detail (all cohorts: S1-S2 mean diff on native log scale)")
    piv = rs.pivot_table(index="gene", columns="cohort", values="diff")
    padj_piv = rs.pivot_table(index="gene", columns="cohort", values="p_adj")
    print("diff (S1-S2)"); print(piv.round(3).to_string())
    print("p_adj")
    print(padj_piv.apply(lambda col: col.map(
        lambda v: f"{v:.2e}" if pd.notna(v) else "NA"
    )).to_string())


if __name__ == "__main__":
    main()