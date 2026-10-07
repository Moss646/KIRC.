"""68-gene NCC classifier cross-cohort concordance, v2.

Uniform mean-difference reference across cohorts. TCGA reference is not on
limma, matching key14_cross_cohort_v2.py and lipid16_cross_cohort_v2.py.

Effect estimator (same for every cohort):
    S1 mean minus S2 mean on the cohort's native log scale.

Cohorts:
    TCGA-KIRC    log2TPM, frozen subtype assignment
    E-MTAB-1980  microarray log intensities, nearest-centroid on z-scored genes
    ICGC RECA-EU log2(normalized_read_count + 1), same assignment
    CPTAC        RNA-seq, same assignment, effect on the raw matrix

Pearson r: external effect vs TCGA mean difference. Fisher-z 95% CI.
Sign convention: >0 means higher in S1.

Output:
    results/tables/classifier68_genes_comparison_v2.csv
    results/tables/classifier68_cross_cohort_summary_v2.csv
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
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from config import (
    GEO_EMTAB_EXPR,
    GEO_NM2GENE,
    GENE_MAPPING_JSON,
    ICGC_EXP_SEQ,
    ICGC_SPECIMEN,
    CPTAC_RNA,
    SUBTYPE_CLASSIFIER_JSON,
    TCGA_EXPR,
    SUBTYPE_ASSIGNMENT_CSV,
    TABLES_DIR,
    emit_provenance,
)


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


def fc_map_cols(mat, s1i, s2i, genes):
    out = {}
    for g in genes:
        if g not in mat.columns:
            continue
        v1 = mat.iloc[s1i][g].astype(float).dropna()
        v2 = mat.iloc[s2i][g].astype(float).dropna()
        if len(v1) >= 3 and len(v2) >= 3:
            out[g] = v1.mean() - v2.mean()
    return out


def tcga_mean_diff(ref_genes):
    expr = pd.read_csv(TCGA_EXPR, index_col=0)
    sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)

    s1 = [c for c in sub[sub.Subtype == "Subtype_1"]["Patient"] if c in expr.columns]
    s2 = [c for c in sub[sub.Subtype == "Subtype_2"]["Patient"] if c in expr.columns]

    return {
        g: float(expr.loc[g, s1].astype(float).mean() - expr.loc[g, s2].astype(float).mean())
        for g in ref_genes
        if g in expr.index
    }


def main():
    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    ref_genes = clf["genes"]
    print(f"Classifier genes: {len(ref_genes)}")

    tcga_fc = tcga_mean_diff(ref_genes)
    print(f"TCGA (mean diff): {len(tcga_fc)}/68 genes")

    emat = prepare_emtab()
    e_s1i, e_s2i = assign_subtypes(emat, ref_genes, clf)
    print(f"E-MTAB: {emat.shape[0]} samples, S1={len(e_s1i)} S2={len(e_s2i)}")
    em_fc = fc_map_cols(emat, e_s1i, e_s2i, ref_genes)
    print(f"  effects computed: {len(em_fc)}/68")

    g2e = load_mapping()
    iem = prepare_icgc(g2e, list(set(ref_genes)))
    i_s1i, i_s2i = assign_subtypes(iem, ref_genes, clf)
    print(f"ICGC:   {iem.shape[0]} samples, S1={len(i_s1i)} S2={len(i_s2i)}")
    ic_fc = fc_map_cols(iem, i_s1i, i_s2i, ref_genes)
    print(f"  effects computed: {len(ic_fc)}/68")

    rnaseq_cp = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
    rnaz = rnaseq_cp.subtract(rnaseq_cp.mean(axis=1), axis=0).div(
        rnaseq_cp.std(axis=1, ddof=1), axis=0
    ).fillna(0)

    cp_s1i, cp_s2i = assign_subtypes(rnaz.T, ref_genes, clf)
    print(f"CPTAC:  {rnaseq_cp.shape[1]} samples, S1={len(cp_s1i)} S2={len(cp_s2i)}")

    cp_s1_cols = rnaseq_cp.columns[cp_s1i]
    cp_s2_cols = rnaseq_cp.columns[cp_s2i]

    cp_fc = {}
    for g in ref_genes:
        if g not in rnaseq_cp.index:
            continue
        v1 = rnaseq_cp.loc[g, cp_s1_cols].astype(float).dropna()
        v2 = rnaseq_cp.loc[g, cp_s2_cols].astype(float).dropna()
        if len(v1) >= 3 and len(v2) >= 3:
            cp_fc[g] = v1.mean() - v2.mean()
    print(f"  effects computed: {len(cp_fc)}/68")

    print("\n" + "=" * 60)
    print("PEARSON r vs TCGA mean difference")
    print("=" * 60)

    ref_sign = {g: (1 if tcga_fc[g] > 0 else -1) for g in tcga_fc}
    summary_rows = []
    comp_rows = []

    for g, v in sorted(tcga_fc.items()):
        comp_rows.append({
            "gene": g,
            "cohort": "TCGA",
            "diff": v,
            "effect_measure": "mean_diff_log2TPM",
        })

    for name, fc in [("E-MTAB", em_fc), ("ICGC", ic_fc), ("CPTAC", cp_fc)]:
        common = sorted(set(tcga_fc) & set(fc))
        x = [tcga_fc[g] for g in common]
        y = [fc[g] for g in common]

        r_val, p_val = pearsonr(x, y)
        z = np.arctanh(r_val)
        se = 1 / np.sqrt(len(x) - 3)
        r_lo = np.tanh(z - 1.96 * se)
        r_hi = np.tanh(z + 1.96 * se)

        dir_match = sum(
            1 for g in common if (1 if fc[g] > 0 else -1) == ref_sign[g]
        )

        print(
            f"  {name:6s}: n={len(common):2d}  r={r_val:+.3f} "
            f"({r_lo:.3f}-{r_hi:.3f})  P={p_val:.1e}  "
            f"dir_match={dir_match}/{len(common)}"
        )

        summary_rows.append({
            "cohort": name,
            "n_detected": len(common),
            "dir_match": dir_match,
            "pearson_n": len(common),
            "pearson_r": r_val,
            "pearson_CI": f"{r_lo:.3f}-{r_hi:.3f}",
            "pearson_P": p_val,
        })

        for g in common:
            comp_rows.append({
                "gene": g,
                "cohort": name,
                "diff": fc[g],
                "effect_measure": "mean_diff_log_scale",
            })

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(TABLES_DIR, "classifier68_cross_cohort_summary_v2.csv"),
        index=False,
    )
    pd.DataFrame(comp_rows).to_csv(
        os.path.join(TABLES_DIR, "classifier68_genes_comparison_v2.csv"),
        index=False,
    )
    print("\nSaved: classifier68_cross_cohort_summary_v2.csv, "
          "classifier68_genes_comparison_v2.csv")

    cname = {"E-MTAB": "E-MTAB-1980", "ICGC": "ICGC RECA-EU", "CPTAC": "CPTAC"}
    prov = []
    for s in summary_rows:
        if s["cohort"] not in cname:
            continue
        c = cname[s["cohort"]]
        prov.append({
            "section": "Classifier robustness",
            "item": "68-gene Pearson r vs TCGA",
            "cohort": c,
            "value": f"+{s['pearson_r']:.3f}",
        })
        prov.append({
            "section": "Classifier robustness",
            "item": "68-gene 95% CI",
            "cohort": c,
            "value": s["pearson_CI"],
        })
        prov.append({
            "section": "Classifier robustness",
            "item": "68-gene Pearson P",
            "cohort": c,
            "value": f"{s['pearson_P']:.1e}",
        })

    emit_provenance("classifier_concordance_v2.py", prov)


if __name__ == "__main__":
    main()