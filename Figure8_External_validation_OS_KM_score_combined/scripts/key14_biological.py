"""14-gene lipid-metabolism panel: cross-cohort biological validation.

Block 1: direction concordance and BH-significant gene counts.
Block 2: effect-size Pearson correlation against TCGA-KIRC.

Each block emits provenance under its own script name so the summary
collector contract stays unchanged.

Outputs:
    results/tables/key14_genes_comparison.csv
    results/tables/_provenance/key14_compare.json
    results/tables/_provenance/key14_pearson_all.json
"""

import gzip
import io
import json
import os
import sys
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, pearsonr
from statsmodels.stats.multitest import multipletests

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

from config import *

plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"

C1 = "#C65A5A"
C2 = "#4C78A8"

KEY_GENES = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36",
    "SLC27A2", "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]


def run_compare():
    key_genes = list(KEY_GENES)

    with open(BALANCED_GENES_TXT) as f:
        ncc = [line.strip() for line in f]

    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    tcga_s1 = np.array(clf["centroid_s1"])
    tcga_s2 = np.array(clf["centroid_s2"])

    # E-MTAB-1980
    expr = pd.read_csv(GEO_EMTAB_EXPR, sep="\t", index_col=0)
    sys_col = expr["SystematicName"].values

    with open(GEO_NM2GENE) as f:
        nm2g = json.load(f)

    gene_data = defaultdict(list)
    for i in range(len(sys_col)):
        nm = sys_col[i]
        if nm.startswith("NM_") and nm in nm2g:
            gene_data[nm2g[nm]].append(expr.iloc[i, 2:].astype(float).values)

    emat_em = pd.DataFrame(
        {g: np.mean(arrs, axis=0) for g, arrs in gene_data.items()},
        index=[c.split("_")[0] for c in expr.columns[2:]],
    )
    ez = emat_em.subtract(emat_em.mean(axis=0), axis=1).div(
        emat_em.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    com = [g for g in ncc if g in ez.columns]
    tci = [ncc.index(g) for g in com]
    d1 = cdist(ez[com].values, tcga_s1[tci].reshape(1, -1))[:, 0]
    d2 = cdist(ez[com].values, tcga_s2[tci].reshape(1, -1))[:, 0]

    em_s1 = [emat_em.index[i] for i in range(len(emat_em)) if d1[i] < d2[i]]
    em_s2 = [emat_em.index[i] for i in range(len(emat_em)) if d1[i] >= d2[i]]

    # ICGC RECA-EU
    all_symbols = list(set(ncc + key_genes))
    g2e = get_gene_ensembl_mapping()
    e2s = {v: k for k, v in g2e.items()}
    print(f"Mapped: {len(g2e)}/{len(all_symbols)}")

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
        {d: {g: iexp[d].get(g, 0) for g in all_symbols if g in g2e} for d in iexp}
    ).T.astype(float)
    iem = np.log2(iem + 1)
    iez = iem.subtract(iem.mean(axis=0), axis=1).div(
        iem.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    ico = [g for g in ncc if g in iez.columns]
    iti = [ncc.index(g) for g in ico]
    id1 = cdist(iez[ico].values, tcga_s1[iti].reshape(1, -1))[:, 0]
    id2 = cdist(iez[ico].values, tcga_s2[iti].reshape(1, -1))[:, 0]

    ic_s1 = [iem.index[i] for i in range(len(iem)) if id1[i] < id2[i]]
    ic_s2 = [iem.index[i] for i in range(len(iem)) if id1[i] >= id2[i]]

    all_rs = []
    for g in key_genes:
        if g in emat_em.columns:
            v1 = emat_em.loc[em_s1, g].astype(float)
            v2 = emat_em.loc[em_s2, g].astype(float)
            all_rs.append({
                "cohort": "E-MTAB",
                "gene": g,
                "S1": v1.mean(),
                "S2": v2.mean(),
                "p": mannwhitneyu(v1, v2)[1],
            })
        else:
            all_rs.append({
                "cohort": "E-MTAB",
                "gene": g,
                "S1": np.nan,
                "S2": np.nan,
                "p": np.nan,
            })

        if g in iem.columns:
            v1 = iem.loc[ic_s1, g].astype(float)
            v2 = iem.loc[ic_s2, g].astype(float)
            all_rs.append({
                "cohort": "ICGC",
                "gene": g,
                "S1": v1.mean(),
                "S2": v2.mean(),
                "p": mannwhitneyu(v1, v2)[1],
            })
        else:
            all_rs.append({
                "cohort": "ICGC",
                "gene": g,
                "S1": np.nan,
                "S2": np.nan,
                "p": np.nan,
            })

    tcga_expr = pd.read_csv(TCGA_EXPR, index_col=0)
    tcga_sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)
    tcga_s1p = tcga_sub[tcga_sub["Subtype"] == "Subtype_1"]["Patient"].tolist()
    tcga_s2p = tcga_sub[tcga_sub["Subtype"] == "Subtype_2"]["Patient"].tolist()

    for g in key_genes:
        if g in tcga_expr.index:
            v1 = tcga_expr.loc[g, [c for c in tcga_s1p if c in tcga_expr.columns]].astype(float)
            v2 = tcga_expr.loc[g, [c for c in tcga_s2p if c in tcga_expr.columns]].astype(float)
            all_rs.append({
                "cohort": "TCGA",
                "gene": g,
                "S1": v1.mean(),
                "S2": v2.mean(),
                "p": mannwhitneyu(v1, v2)[1],
            })

    rs = pd.DataFrame(all_rs)
    rs["diff"] = rs["S1"] - rs["S2"]
    rs["logp"] = -np.log10(rs["p"].clip(lower=1e-50))

    rs["p_adj"] = np.nan
    for c in rs["cohort"].unique():
        mask = rs["cohort"] == c
        pvals = rs.loc[mask, "p"].dropna()
        if len(pvals) > 0:
            _, padj, _, _ = multipletests(pvals, method="fdr_bh")
            rs.loc[pvals.index, "p_adj"] = padj

    print("Summary (BH corrected per cohort)")
    for c in ["TCGA", "E-MTAB", "ICGC"]:
        sub = rs[rs.cohort == c]
        n_sig_raw = sum(sub.p < 0.05)
        n_sig_bh = sum(sub.p_adj < 0.05)
        n_total = sub.dropna(subset=["p"]).shape[0]

        if n_sig_bh > 0:
            sig_genes = sub[sub.p_adj < 0.05].gene
            tcga_diff = rs[(rs.cohort == "TCGA") & (rs.gene.isin(sig_genes))]["diff"]
            same_dir = sum(
                (sub[sub.p_adj < 0.05]["diff"] > 0) == (tcga_diff > 0).values
            )
        else:
            same_dir = 0

        print(
            f"{c}: {n_total} genes, raw sig={n_sig_raw}, BH sig={n_sig_bh}, "
            f"same direction vs TCGA in BH-sig genes: {same_dir}/{n_sig_bh}"
        )

    rs.to_csv(os.path.join(TABLES_DIR, "key14_genes_comparison.csv"), index=False)
    print("\nSaved: key14_genes_comparison.csv")

    try:
        tcga_diff = rs[rs.cohort == "TCGA"].set_index("gene")["diff"]
        rows = []

        for c in ["TCGA", "E-MTAB", "ICGC"]:
            sub = rs[rs.cohort == c].dropna(subset=["diff"])
            total = len(sub)

            if c == "TCGA":
                match = total
            else:
                match = int(sum(
                    (sub["diff"] > 0).values == (tcga_diff.reindex(sub["gene"]) > 0).values
                ))

            if "p_adj" in sub:
                bh = int((sub["p_adj"] < 0.05).sum())
            else:
                bh = int((sub["p"] < 0.05).sum())

            cohort_name = {
                "TCGA": "TCGA-KIRC (Training)",
                "E-MTAB": "E-MTAB-1980",
                "ICGC": "ICGC RECA-EU",
            }[c]

            rows.append({
                "section": "Biological validation",
                "item": "14-gene direction concordance",
                "cohort": cohort_name,
                "value": f"{match} of {total}",
            })
            rows.append({
                "section": "Biological validation",
                "item": "14-gene BH-sig genes (MWU FDR<0.05)",
                "cohort": cohort_name,
                "value": f"{bh} of {total}",
            })

        emit_provenance("key14_compare.py", rows)
    except Exception as e:
        print(f"[PROVENANCE] key14_compare emission skipped: {e}")


def run_pearson():
    key14 = list(KEY_GENES)

    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    ref_genes = clf["genes"]
    c1 = np.array(clf["centroid_s1"])
    c2 = np.array(clf["centroid_s2"])

    # TCGA effect size
    tcga_expr = pd.read_csv(TCGA_EXPR, index_col=0)
    tcga_sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)
    t1 = [c for c in tcga_sub[tcga_sub.Subtype == "Subtype_1"]["Patient"] if c in tcga_expr.columns]
    t2 = [c for c in tcga_sub[tcga_sub.Subtype == "Subtype_2"]["Patient"] if c in tcga_expr.columns]

    tcga_fc = {
        g: tcga_expr.loc[g, t1].mean() - tcga_expr.loc[g, t2].mean()
        for g in key14
        if g in tcga_expr.index
    }
    print(f"TCGA: {len(tcga_fc)}/14 genes available")

    # E-MTAB-1980
    expr_em = pd.read_csv(GEO_EMTAB_EXPR, sep="\t", index_col=0)
    sys_col = expr_em["SystematicName"].values

    with open(GEO_NM2GENE) as f:
        nm2g = json.load(f)

    gd = defaultdict(list)
    for i in range(len(sys_col)):
        nm = sys_col[i]
        if nm.startswith("NM_") and nm in nm2g:
            gd[nm2g[nm]].append(expr_em.iloc[i, 2:].astype(float).values)

    emat = pd.DataFrame(
        {g: np.mean(arrs, axis=0) for g, arrs in gd.items()},
        index=[c.split("_")[0] for c in expr_em.columns[2:]],
    )
    ez = emat.subtract(emat.mean(axis=0), axis=1).div(
        emat.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    com_em = [g for g in ref_genes if g in ez.columns]
    tci_em = [ref_genes.index(g) for g in com_em]
    ed1 = cdist(ez[com_em].values, c1[tci_em].reshape(1, -1))[:, 0]
    ed2 = cdist(ez[com_em].values, c2[tci_em].reshape(1, -1))[:, 0]
    e_s1i = np.where(ed1 < ed2)[0]
    e_s2i = np.where(ed1 >= ed2)[0]

    em_fc14 = {
        g: emat.iloc[e_s1i][g].mean() - emat.iloc[e_s2i][g].mean()
        for g in key14
        if g in emat.columns
    }
    cm_em = sorted(set(tcga_fc) & set(em_fc14))
    x_em = [tcga_fc[g] for g in cm_em]
    y_em = [em_fc14[g] for g in cm_em]
    r_em, p_em = pearsonr(x_em, y_em)
    z_em = np.arctanh(r_em)
    se_em = 1 / np.sqrt(len(x_em) - 3)
    print(
        f"E-MTAB: n={len(cm_em)} r={r_em:.3f} "
        f"({np.tanh(z_em - 1.96 * se_em):.3f}\u2013{np.tanh(z_em + 1.96 * se_em):.3f}) "
        f"P={p_em:.1e}"
    )

    # ICGC RECA-EU
    all_sym = list(set(key14 + list(ref_genes)))
    g2e = get_gene_ensembl_mapping()
    e2s = {v: k for k, v in g2e.items()}

    tds = set()
    with gzip.open(ICGC_SPECIMEN, "rt") as f:
        hdr = f.readline().strip().split("\t")
        sct = hdr.index("specimen_type")
        did = hdr.index("icgc_donor_id")
        for line in f:
            c = line.strip().split("\t")
            if "Primary tumour" in c[sct]:
                tds.add(c[did])

    iexp = defaultdict(lambda: defaultdict(float))
    with gzip.open(ICGC_EXP_SEQ, "rt") as f:
        hdr = f.readline().strip().split("\t")
        gi = hdr.index("gene_id")
        si = hdr.index("icgc_donor_id")
        fi = hdr.index("normalized_read_count")
        for line in f:
            c = line.strip().split("\t")
            gv = c[gi].replace("gene:", "")
            if gv in e2s and c[si] in tds:
                iexp[c[si]][e2s[gv]] = float(c[fi])

    iem = pd.DataFrame(
        {d: {g: iexp[d].get(g, 0) for g in all_sym if g in g2e} for d in iexp}
    ).T.astype(float)
    iem = np.log2(iem + 1)
    iez = iem.subtract(iem.mean(axis=0), axis=1).div(
        iem.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    com_ic = [g for g in ref_genes if g in iez.columns]
    tci_ic = [ref_genes.index(g) for g in com_ic]
    id1 = cdist(iez[com_ic].values, c1[tci_ic].reshape(1, -1))[:, 0]
    id2 = cdist(iez[com_ic].values, c2[tci_ic].reshape(1, -1))[:, 0]
    i_s1i = np.where(id1 < id2)[0]
    i_s2i = np.where(id1 >= id2)[0]

    ic_fc14 = {
        g: iem.iloc[i_s1i][g].mean() - iem.iloc[i_s2i][g].mean()
        for g in key14
        if g in iem.columns
    }
    cm_ic = sorted(set(tcga_fc) & set(ic_fc14))
    x_ic = [tcga_fc[g] for g in cm_ic]
    y_ic = [ic_fc14[g] for g in cm_ic]
    r_ic, p_ic = pearsonr(x_ic, y_ic)
    z_ic = np.arctanh(r_ic)
    se_ic = 1 / np.sqrt(len(x_ic) - 3)
    print(
        f"ICGC:  n={len(cm_ic)} r={r_ic:.3f} "
        f"({np.tanh(z_ic - 1.96 * se_ic):.3f}\u2013{np.tanh(z_ic + 1.96 * se_ic):.3f}) "
        f"P={p_ic:.1e}"
    )

    # CPTAC protein
    prot = pd.read_csv(CPTAC_PROTEIN, sep="\t", index_col=0)
    rnaseq_cp = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
    rnaz_cp = rnaseq_cp.subtract(rnaseq_cp.mean(axis=1), axis=0).div(
        rnaseq_cp.std(axis=1, ddof=1), axis=0
    ).fillna(0)

    com_cp = [g for g in ref_genes if g in rnaz_cp.index]
    tci_cp = [ref_genes.index(g) for g in com_cp]
    cd1 = cdist(rnaz_cp.loc[com_cp].T.values, c1[tci_cp].reshape(1, -1))[:, 0]
    cd2 = cdist(rnaz_cp.loc[com_cp].T.values, c2[tci_cp].reshape(1, -1))[:, 0]
    cp_s1s = rnaseq_cp.columns[cd1 < cd2]
    cp_s2s = rnaseq_cp.columns[cd1 >= cd2]

    cp_prot_fc = {}
    for g in key14:
        if g not in prot.index:
            continue
        v1 = prot.loc[g, [c for c in cp_s1s if c in prot.columns]].astype(float)
        v2 = prot.loc[g, [c for c in cp_s2s if c in prot.columns]].astype(float)
        v1 = v1[~np.isnan(v1)]
        v2 = v2[~np.isnan(v2)]
        if len(v1) >= 3 and len(v2) >= 3:
            cp_prot_fc[g] = v1.mean() - v2.mean()

    cm_cp = sorted(set(tcga_fc) & set(cp_prot_fc))
    x_cp = [tcga_fc[g] for g in cm_cp]
    y_cp = [cp_prot_fc[g] for g in cm_cp]
    r_cp, p_cp = pearsonr(x_cp, y_cp)
    z_cp = np.arctanh(r_cp)
    se_cp = 1 / np.sqrt(len(x_cp) - 3)
    print(
        f"CPTAC: n={len(cm_cp)} r={r_cp:.3f} "
        f"({np.tanh(z_cp - 1.96 * se_cp):.3f}\u2013{np.tanh(z_cp + 1.96 * se_cp):.3f}) "
        f"P={p_cp:.1e}"
    )

    print("\n" + "=" * 60)
    print("TABLE COMPARISON")
    print("=" * 60)

    for name, r_val, r_tab, p_val, p_tab in [
        ("E-MTAB", r_em, 0.940, p_em, 5.5e-6),
        ("ICGC", r_ic, 0.740, p_ic, 2.5e-3),
        ("CPTAC", r_cp, 0.848, p_cp, 4.9e-4),
    ]:
        ok_r = "OK" if abs(r_val - r_tab) < 0.02 else f"MISMATCH ({r_val:.3f})"
        ok_p = "OK" if abs(p_val - p_tab) / p_tab < 0.1 else f"MISMATCH ({p_val:.1e})"
        print(f"  {name}: r={ok_r} P={ok_p}")

    try:
        def ci_str(rv, zv, sev):
            return f"{np.tanh(zv - 1.96 * sev):.3f}\u2013{np.tanh(zv + 1.96 * sev):.3f}"

        rows = []
        cohorts = [
            ("E-MTAB-1980", r_em, p_em, z_em, se_em),
            ("ICGC RECA-EU", r_ic, p_ic, z_ic, se_ic),
            ("CPTAC", r_cp, p_cp, z_cp, se_cp),
        ]
        for name, rv, pv, zv, sev in cohorts:
            rows.append({
                "section": "Biological validation",
                "item": "14-gene Pearson r vs TCGA",
                "cohort": name,
                "value": f"+{rv:.3f}",
            })
            rows.append({
                "section": "Biological validation",
                "item": "14-gene 95% CI",
                "cohort": name,
                "value": ci_str(rv, zv, sev),
            })
            rows.append({
                "section": "Biological validation",
                "item": "14-gene Pearson P",
                "cohort": name,
                "value": f"{pv:.1e}",
            })

        emit_provenance("key14_pearson_all.py", rows)
    except Exception as e:
        print(f"[PROVENANCE] key14_pearson_all emission skipped: {e}")


def main():
    run_compare()
    run_pearson()


if __name__ == "__main__":
    main()