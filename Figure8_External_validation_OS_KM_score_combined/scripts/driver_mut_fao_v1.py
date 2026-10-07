"""Driver mutation vs FAO-low analysis.

Outputs:
    results/tables/driver_mut_fao_v1.csv
    results/tables/driver_gistic_del_v1.csv
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from config import (
    DATA_DIR,
    TABLES_DIR,
    SUBTYPE_ASSIGNMENT_CSV,
)

GENES = ["VHL", "PBRM1", "BAP1", "SETD2", "KDM5C", "MTOR", "TP53"]

NONSYN = [
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "In_Frame_Del",
    "In_Frame_Ins",
]


def _find(name):
    """Look for a file under DATA_DIR (recursive)."""
    direct = DATA_DIR / name
    if direct.exists():
        return direct
    hits = list(Path(DATA_DIR).rglob(name))
    if hits:
        return hits[0]
    return direct


MAF = _find("TCGA-KIRC.merged.maf")
GISTIC = _find("TCGA-KIRC.gistic.tsv")
SUBTYPE_CSV = Path(SUBTYPE_ASSIGNMENT_CSV)

FAO_SCORE_CSV = Path(TABLES_DIR) / "nonncc9_fao_score_per_sample_v1.csv"


def bh_correct(pvals):
    arr = np.asarray(pvals, dtype=float)
    n = len(arr)
    order = np.argsort(arr)
    sorted_p = arr[order]
    adj = sorted_p * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out = np.zeros(n)
    out[order] = adj
    return out


def main():
    for p in [MAF, GISTIC, SUBTYPE_CSV, FAO_SCORE_CSV]:
        if not Path(p).exists():
            raise SystemExit(f"missing input: {p}")

    # Subtype map
    sub = pd.read_csv(SUBTYPE_CSV)
    sub["S"] = sub["Subtype"].map({"Subtype_1": "S1", "Subtype_2": "S2"})
    sub_map = dict(zip(sub["Patient"], sub["S"]))
    n_s1 = sum(v == "S1" for v in sub_map.values())
    n_s2 = sum(v == "S2" for v in sub_map.values())
    print(f"subtype map: S1={n_s1}, S2={n_s2}")

    # MAF
    maf = pd.read_csv(MAF, sep="\t", comment="#", low_memory=False)
    maf = maf[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"]].copy()
    maf = maf[maf["Variant_Classification"].isin(NONSYN)].copy()
    maf["Patient"] = maf["Tumor_Sample_Barcode"].str[:12]

    mut_pats = sorted(set(maf["Patient"]) & set(sub_map))
    s1_pats = [p for p in mut_pats if sub_map[p] == "S1"]
    s2_pats = [p for p in mut_pats if sub_map[p] == "S2"]
    print(
        f"MAF nonsyn variants: {len(maf)}; mutation cohort: "
        f"{len(mut_pats)} (S1={len(s1_pats)}, S2={len(s2_pats)})"
    )

    # FAO9 per-sample
    fao9_all = pd.read_csv(FAO_SCORE_CSV)
    fao9 = fao9_all[fao9_all["cohort"] == "TCGA"][["sample", "score"]].copy()
    fao9 = fao9.rename(columns={"score": "FAO9"})
    fao9["Patient"] = fao9["sample"].str[:12]
    fao9_map = dict(zip(fao9["Patient"], fao9["FAO9"]))
    print(f"FAO9 scores available: {len(fao9_map)}")

    # A + B
    rows = []
    for gene in GENES:
        g = maf[maf["Hugo_Symbol"] == gene]
        s1_mut = set(g["Patient"]) & set(s1_pats)
        s2_mut = set(g["Patient"]) & set(s2_pats)
        n1, n2 = len(s1_pats), len(s2_pats)
        m1, m2 = len(s1_mut), len(s2_mut)
        orr, p_fish = fisher_exact([[m1, n1 - m1], [m2, n2 - m2]])

        mut_set = set(g["Patient"])
        wt_pats = [p for p in mut_pats if p not in mut_set]
        v_mut = [fao9_map[p] for p in mut_pats if p in mut_set and p in fao9_map]
        v_wt = [fao9_map[p] for p in wt_pats if p in fao9_map]

        if len(v_mut) >= 3 and len(v_wt) >= 3:
            u, p_mwu = mannwhitneyu(v_mut, v_wt, alternative="two-sided")
            med_mut = pd.Series(v_mut).median()
            med_wt = pd.Series(v_wt).median()
        else:
            p_mwu = float("nan")
            med_mut = float("nan")
            med_wt = float("nan")

        rows.append(dict(
            gene=gene,
            S1_mut=m1, S1_n=n1, S1_pct=100 * m1 / n1,
            S2_mut=m2, S2_n=n2, S2_pct=100 * m2 / n2,
            OR_S1vsS2=orr, p_fisher=p_fish,
            FAO9_n_mut=len(v_mut), FAO9_n_wt=len(v_wt),
            FAO9_med_mut=med_mut, FAO9_med_wt=med_wt,
            p_mwu_FAO9=p_mwu,
        ))

    res = pd.DataFrame(rows)
    res["p_bh_fisher"] = bh_correct(res["p_fisher"].tolist())
    res["p_bh_mwu"] = bh_correct(res["p_mwu_FAO9"].tolist())
    res.to_csv(Path(TABLES_DIR) / "driver_mut_fao_v1.csv", index=False)
    print("\nA+B: driver mutations vs S1/S2 + FAO9 score")
    print(res.round(4).to_string(index=False))

    # C: GISTIC
    gist = pd.read_csv(GISTIC, sep="\t", index_col=0)
    rows_c = []
    for gene in ["VHL", "PBRM1", "BAP1", "SETD2"]:
        if gene not in gist.index:
            print(f"  [warn] {gene} not in GISTIC table")
            continue

        row = gist.loc[gene]
        samples = gist.columns
        pat = pd.Series(samples, index=samples).str[:12]

        df = pd.DataFrame({"Patient": pat.values, "val": row.values})
        df = df.groupby("Patient")["val"].min().rename("val").reset_index()
        df = df[df["Patient"].isin(sub_map)].copy()
        df["S"] = df["Patient"].map(sub_map)

        n_all = len(df)
        del_all = int((df["val"] <= -1).sum())
        s1 = df[df["S"] == "S1"]
        s2 = df[df["S"] == "S2"]
        d1, d2 = int((s1["val"] <= -1).sum()), int((s2["val"] <= -1).sum())
        _, p_del = fisher_exact([[d1, len(s1) - d1], [d2, len(s2) - d2]])

        rows_c.append(dict(
            gene=gene,
            n=n_all,
            del_n=del_all, del_pct=100 * del_all / n_all,
            S1_del=d1, S1_n=len(s1), S1_del_pct=100 * d1 / len(s1),
            S2_del=d2, S2_n=len(s2), S2_del_pct=100 * d2 / len(s2),
            p_fisher_del_S1vsS2=p_del,
        ))

    res_c = pd.DataFrame(rows_c)
    res_c.to_csv(Path(TABLES_DIR) / "driver_gistic_del_v1.csv", index=False)
    print("\nC: GISTIC gene-level deletion rates (3p drivers)")
    print(res_c.round(4).to_string(index=False))


if __name__ == "__main__":
    main()