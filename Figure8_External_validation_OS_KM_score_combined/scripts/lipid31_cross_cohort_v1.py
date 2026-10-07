"""Combined 30-gene cross-cohort analysis: key14 panel + lipid remodeling 16.

The two panels are disjoint (0 overlap):
    key14   FAO core, uptake, synthesis, storage
    lipid16 phospholipid remodeling enzymes (16 genes, no PLIN2, no PLAAT3)

Results are reported per panel and combined.

Effect estimator (same for every cohort):
    S1 mean minus S2 mean on the cohort's native log scale.

Cohorts:
    TCGA-KIRC    log2TPM, frozen subtype assignment, MWU + BH
    E-MTAB-1980  microarray log intensities, per-gene MWU + BH
    ICGC RECA-EU log2(normalized_read_count + 1), per-gene MWU + BH
    CPTAC        protein log abundance, per-gene MWU + BH

Subtyping: frozen 68-gene nearest-centroid NCC in every cohort.
Sign convention: >0 means higher in S1.

No provenance is emitted. Tables 1-3 build chain is untouched.

Outputs:
    results/tables/lipid31_genes_comparison_v1.csv
    results/tables/lipid31_cross_cohort_summary_v1.csv
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
    SUBTYPE_CLASSIFIER_JSON,
    BALANCED_GENES_TXT,
    GEO_EMTAB_EXPR,
    GEO_NM2GENE,
    GENE_MAPPING_JSON,
    ICGC_EXP_SEQ,
    ICGC_SPECIMEN,
    CPTAC_PROTEIN,
    CPTAC_RNA,
    TCGA_EXPR,
    SUBTYPE_ASSIGNMENT_CSV,
    TABLES_DIR,
)

KEY14 = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36", "SLC27A2",
    "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]
LIPID16 = [
    "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "MBOAT1", "MBOAT2", "LCLAT1", "MBOAT7",
    "PLA2G4A", "PLA2G6", "PNPLA8", "PLA2G4C",
    "PLA2G2A", "PLA2G4F", "PLA2G1B", "PLA2G2D",
]
ALL31 = KEY14 + LIPID16
assert len(ALL31) == 30 and len(set(ALL31)) == 30, "gene sets must be disjoint"

GROUP = {g: "key14_panel" for g in KEY14}
for g in LIPID16:
    GROUP[g] = "lipid16_tumor_intrinsic"

LIPID16_ENSEMBL_JSON = os.path.join(
    os.path.dirname(SCRIPT_DIR), "data", "lipid16_ensembl_mapping_v1.json"
)


def load_centroids():
    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)
    with open(BALANCED_GENES_TXT) as f:
        ncc = [line.strip() for line in f if line.strip()]
    return clf, ncc


def nearest_centroid(zmat, clf, ncc):
    c1 = np.array(clf["centroid_s1"])
    c2 = np.array(clf["centroid_s2"])
    com = [g for g in ncc if g in zmat.columns]
    tci = [ncc.index(g) for g in com]
    d1 = cdist(zmat[com].values, c1[tci].reshape(1, -1))[:, 0]
    d2 = cdist(zmat[com].values, c2[tci].reshape(1, -1))[:, 0]
    return list(np.where(d1 < d2)[0]), list(np.where(d1 >= d2)[0])


def prepare_emtab(clf, ncc):
    expr = pd.read_csv(GEO_EMTAB_EXPR, sep="\t", index_col=0)
    sys_col = expr["SystematicName"].values

    with open(GEO_NM2GENE) as f:
        nm2g = json.load(f)

    gene_data = defaultdict(list)
    for i in range(len(sys_col)):
        nm = sys_col[i]
        if nm.startswith("NM_") and nm in nm2g:
            gene_data[nm2g[nm]].append(expr.iloc[i, 2:].astype(float).values)

    emat = pd.DataFrame(
        {g: np.mean(arrs, axis=0) for g, arrs in gene_data.items()},
        index=[c.split("_")[0] for c in expr.columns[2:]],
    )
    ez = emat.subtract(emat.mean(axis=0), axis=1).div(
        emat.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    s1i, s2i = nearest_centroid(ez, clf, ncc)
    return emat, s1i, s2i


def prepare_icgc(clf, ncc, need_genes):
    with open(LIPID16_ENSEMBL_JSON) as f:
        lipid_map = json.load(f)
    with open(GENE_MAPPING_JSON) as f:
        panel_map = json.load(f)

    g2e = dict(panel_map)
    g2e.update(lipid_map)
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

    need = set(ncc) | set(need_genes)
    iem = pd.DataFrame(
        {d: {g: iexp[d].get(g, 0) for g in need if g in g2e} for d in iexp}
    ).T.astype(float)
    iem = np.log2(iem + 1)
    iez = iem.subtract(iem.mean(axis=0), axis=1).div(
        iem.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    s1i, s2i = nearest_centroid(iez, clf, ncc)
    return iem, s1i, s2i


def external_rows(cohort, mat, s1i, s2i, genes):
    rows = []
    for g in genes:
        if g in mat.columns:
            v1 = mat.iloc[s1i][g].astype(float).dropna()
            v2 = mat.iloc[s2i][g].astype(float).dropna()
            if len(v1) >= 3 and len(v2) >= 3:
                rows.append({
                    "cohort": cohort,
                    "gene": g,
                    "group": GROUP[g],
                    "S1": v1.mean(),
                    "S2": v2.mean(),
                    "p": mannwhitneyu(v1, v2)[1],
                    "effect_measure": "mean_diff_log_scale",
                })
                continue

        rows.append({
            "cohort": cohort,
            "gene": g,
            "group": GROUP[g],
            "S1": np.nan,
            "S2": np.nan,
            "p": np.nan,
            "effect_measure": "mean_diff_log_scale",
        })
    return rows


def tcga_mean_diff_rows():
    expr = pd.read_csv(TCGA_EXPR, index_col=0)
    sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)

    s1 = [c for c in sub[sub.Subtype == "Subtype_1"]["Patient"] if c in expr.columns]
    s2 = [c for c in sub[sub.Subtype == "Subtype_2"]["Patient"] if c in expr.columns]

    rows = []
    for g in ALL31:
        if g in expr.index:
            v1 = expr.loc[g, s1].astype(float)
            v2 = expr.loc[g, s2].astype(float)
            rows.append({
                "cohort": "TCGA",
                "gene": g,
                "group": GROUP[g],
                "S1": float(v1.mean()),
                "S2": float(v2.mean()),
                "p": float(mannwhitneyu(v1, v2)[1]),
                "diff": float(v1.mean() - v2.mean()),
                "effect_measure": "mean_diff_log2TPM",
            })
        else:
            rows.append({
                "cohort": "TCGA",
                "gene": g,
                "group": GROUP[g],
                "S1": np.nan,
                "S2": np.nan,
                "p": np.nan,
                "diff": np.nan,
                "effect_measure": "mean_diff_log2TPM",
            })
    return rows


def cptac_rows(clf, ncc):
    prot = pd.read_csv(CPTAC_PROTEIN, sep="\t", index_col=0)
    rnaseq = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
    rnaz = rnaseq.subtract(rnaseq.mean(axis=1), axis=0).div(
        rnaseq.std(axis=1, ddof=1), axis=0
    ).fillna(0)

    s1i, s2i = nearest_centroid(rnaz.T, clf, ncc)
    s1c = [c for c in rnaseq.columns[s1i] if c in prot.columns]
    s2c = [c for c in rnaseq.columns[s2i] if c in prot.columns]
    print(f"CPTAC:  protein S1={len(s1c)} S2={len(s2c)}")

    rows = []
    for g in ALL31:
        if g in prot.index:
            v1 = prot.loc[g, s1c].astype(float).dropna()
            v2 = prot.loc[g, s2c].astype(float).dropna()
            if len(v1) >= 3 and len(v2) >= 3:
                rows.append({
                    "cohort": "CPTAC",
                    "gene": g,
                    "group": GROUP[g],
                    "S1": v1.mean(),
                    "S2": v2.mean(),
                    "p": mannwhitneyu(v1, v2)[1],
                    "effect_measure": "protein_mean_diff",
                })
                continue

        rows.append({
            "cohort": "CPTAC",
            "gene": g,
            "group": GROUP[g],
            "S1": np.nan,
            "S2": np.nan,
            "p": np.nan,
            "effect_measure": "protein_mean_diff",
        })
    return rows


def main():
    clf, ncc = load_centroids()

    emat, em_s1, em_s2 = prepare_emtab(clf, ncc)
    print(f"E-MTAB: {emat.shape[0]} samples, S1={len(em_s1)} S2={len(em_s2)}")

    iem, ic_s1, ic_s2 = prepare_icgc(clf, ncc, ALL31)
    print(f"ICGC:   {iem.shape[0]} samples, S1={len(ic_s1)} S2={len(ic_s2)}")

    rows = tcga_mean_diff_rows()
    rows += external_rows("E-MTAB", emat, em_s1, em_s2, ALL31)
    rows += external_rows("ICGC", iem, ic_s1, ic_s2, ALL31)
    rows += cptac_rows(clf, ncc)

    rs = pd.DataFrame(rows)
    rs["p_adj"] = rs.get("p_adj", np.nan)
    ext = rs["cohort"] != "TCGA"
    rs.loc[ext, "diff"] = rs.loc[ext, "S1"] - rs.loc[ext, "S2"]

    for c in ["TCGA", "E-MTAB", "ICGC", "CPTAC"]:
        m = rs["cohort"] == c
        pv = rs.loc[m, "p"].dropna()
        if len(pv) > 0:
            _, padj, _, _ = multipletests(pv, method="fdr_bh")
            rs.loc[pv.index, "p_adj"] = padj

    tcga_dir = rs[rs.cohort == "TCGA"].set_index("gene")["diff"]

    print("\nDirection concordance (sign of S1-S2 vs TCGA mean diff)")
    summary = []
    n_total = len(ALL31)

    def concord_block(cohort, panel_name, gene_subset):
        sub = rs[(rs.cohort == cohort) & (rs.gene.isin(gene_subset))].dropna(subset=["diff"])
        if len(sub) == 0:
            return

        match = int(
            ((sub["diff"] > 0).values == (tcga_dir.reindex(sub["gene"]) > 0).values).sum()
        )
        bh = sub[sub.p_adj < 0.05]

        if len(bh):
            bh_match = int(
                ((bh["diff"] > 0).values == (tcga_dir.reindex(bh["gene"]) > 0).values).sum()
            )
        else:
            bh_match = 0

        print(
            f"{cohort:6s} [{panel_name:10s}]: {len(sub)}/{len(gene_subset)} detected, "
            f"dir match {match}/{len(sub)}, BH-sig {len(bh)} (match {bh_match}/{len(bh)})"
        )
        summary.append({
            "cohort": cohort,
            "panel": panel_name,
            "n_detected": len(sub),
            "dir_match": match,
            "n_bh_sig": len(bh),
            "bh_dir_match": bh_match,
        })

    tcga_all = rs[rs.cohort == "TCGA"].dropna(subset=["diff"])
    tcga_sig = tcga_all[tcga_all.p_adj < 0.05]
    print(f"TCGA  : {len(tcga_all)}/{n_total} genes in log2TPM matrix, MWU BH<0.05: {len(tcga_sig)}")

    for pname, subset in [("all31", ALL31), ("key14", KEY14), ("lipid16", LIPID16)]:
        s = tcga_all[tcga_all.gene.isin(subset)]
        print(f"        {pname}: MWU BH<0.05 {len(s[s.p_adj < 0.05])}/{len(s)}")

    summary.append({
        "cohort": "TCGA",
        "panel": "all31",
        "n_detected": len(tcga_all),
        "dir_match": np.nan,
        "n_bh_sig": len(tcga_sig),
        "bh_dir_match": np.nan,
    })

    for c in ["E-MTAB", "ICGC", "CPTAC"]:
        for pname, subset in [("all31", ALL31), ("key14", KEY14), ("lipid16", LIPID16)]:
            concord_block(c, pname, subset)

    print("\nEffect-size Pearson r vs TCGA mean difference")
    for c in ["E-MTAB", "ICGC", "CPTAC"]:
        for pname, subset in [("all31", ALL31), ("key14", KEY14), ("lipid16", LIPID16)]:
            sub = rs[(rs.cohort == c) & (rs.gene.isin(subset))].dropna(subset=["diff"])
            common = sorted(set(sub.gene) & set(tcga_dir.index))
            if len(common) < 5:
                continue

            x = [tcga_dir[g] for g in common]
            y = [sub.set_index("gene").loc[g, "diff"] for g in common]

            r, p = pearsonr(x, y)
            z = np.arctanh(r)
            se = 1 / np.sqrt(len(common) - 3)
            lo = np.tanh(z - 1.96 * se)
            hi = np.tanh(z + 1.96 * se)

            print(f"{c:6s} [{pname:7s}]: n={len(common)} r={r:.3f} ({lo:.3f}-{hi:.3f}) P={p:.1e}")
            summary.append({
                "cohort": c,
                "panel": pname,
                "pearson_n": len(common),
                "pearson_r": r,
                "pearson_CI": f"{lo:.3f}-{hi:.3f}",
                "pearson_P": p,
            })

    out = os.path.join(TABLES_DIR, "lipid31_genes_comparison_v1.csv")
    rs.to_csv(out, index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(TABLES_DIR, "lipid31_cross_cohort_summary_v1.csv"), index=False)
    print(f"\nSaved: {out}")
    print(f"Saved: {os.path.join(TABLES_DIR, 'lipid31_cross_cohort_summary_v1.csv')}")

    print("\nPer-gene detail (all cohorts: S1-S2 mean diff on native log scale)")
    piv = rs.pivot_table(index=["group", "gene"], columns="cohort", values="diff")
    padj_piv = rs.pivot_table(index=["group", "gene"], columns="cohort", values="p_adj")
    print("diff (S1-S2; >0 = higher in S1)")
    print(piv.round(3).to_string())
    print("p_adj (BH within cohort)")
    print(padj_piv.apply(lambda col: col.map(
        lambda v: f"{v:.2e}" if pd.notna(v) else "NA"
    )).to_string())


if __name__ == "__main__":
    main()