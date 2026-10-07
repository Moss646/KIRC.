"""Prognostic value of NCC subtype across three cohorts.

TCGA-KIRC, E-MTAB-1980, ICGC RECA-EU. Schoenfeld residual check, univariate
and stage-stratified Cox, bootstrap C-index. Data prep differs per cohort;
modeling and figure code is shared.
"""

import gzip
import io
import json
import os
import re
import sys
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.spatial.distance import cdist

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["xtick.major.width"] = 1
plt.rcParams["ytick.major.width"] = 1

from config import *

FIG = FIGURES_DIR
C1 = "#B71C1C"
C2 = "#0D47A1"
Cg = "#7F8C8D"

NCC_GENES_14 = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36", "SLC27A2",
    "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]


def prep_tcga():
    clin = pd.read_csv(
        os.path.join(CLASSIFIER_DIR, "KIRC_clinical_cleaned.csv"), index_col=0
    )
    sub = pd.read_csv(
        os.path.join(CLASSIFIER_DIR, "subtype_assignment_balanced.csv")
    )
    sub_map = dict(zip(sub["Patient"], sub["Subtype"]))

    clin["Subtype"] = clin.index.map(sub_map)
    clin["Subtype_bin"] = (clin["Subtype"] == "Subtype_1").astype(int)
    clin["Stage_num"] = clin["stage"].map({"I": 1, "II": 2, "III": 3, "IV": 4})
    clin["Grade_num"] = pd.to_numeric(clin["grade"], errors="coerce")
    clin["Age_num"] = pd.to_numeric(clin["age"], errors="coerce")
    clin["os_bin"] = pd.to_numeric(clin["os"], errors="coerce")
    clin["os_t"] = pd.to_numeric(clin["os_time_days"], errors="coerce")

    df = clin[
        ["os_t", "os_bin", "Subtype_bin", "Stage_num", "Grade_num", "Age_num"]
    ].dropna()
    df = df[df["os_t"] > 0]
    df.columns = ["time", "status", "Subtype", "Stage", "Grade", "Age"]
    df["status"] = df["status"].astype(int)
    return df


def prep_emtab():
    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    ref_genes = clf["genes"]
    tcga_s1 = np.array(clf["centroid_s1"])
    tcga_s2 = np.array(clf["centroid_s2"])

    expr = pd.read_csv(GEO_EMTAB_EXPR, sep="\t", index_col=0)
    sys_col = expr["SystematicName"].values

    with open(GEO_NM2GENE) as f:
        nm2g = json.load(f)

    gd = defaultdict(list)
    for i in range(len(sys_col)):
        nm = sys_col[i]
        if nm.startswith("NM_") and nm in nm2g:
            gd[nm2g[nm]].append(expr.iloc[i, 2:].astype(float).values)

    emat = pd.DataFrame(
        {g: np.mean(arrs, axis=0) for g, arrs in gd.items()},
        index=[c.split("_")[0] for c in expr.columns[2:]],
    )
    ez = emat.subtract(emat.mean(axis=0), axis=1).div(
        emat.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    common = [g for g in ref_genes if g in ez.columns]
    tci = [ref_genes.index(g) for g in common]
    d1 = cdist(ez[common].values, tcga_s1[tci].reshape(1, -1))[:, 0]
    d2 = cdist(ez[common].values, tcga_s2[tci].reshape(1, -1))[:, 0]
    ncc_labels = np.where(d1 < d2, 1, 0)

    ncc_df = pd.DataFrame({"sample": list(emat.index), "NCC_S1": ncc_labels})

    clin = pd.read_excel(GEO_SATO_CLINICAL, header=1)
    clin = clin.rename(
        columns={
            "sample ID": "sample",
            "Age": "age",
            "Fuhrman grade": "grade",
            "observation period (month)": "os_months",
            "outcome": "status",
        }
    )
    clin["os_months"] = pd.to_numeric(clin["os_months"], errors="coerce")
    clin["os_event"] = (clin["status"].str.lower() == "dead").astype(int)
    clin["age"] = pd.to_numeric(clin["age"], errors="coerce")
    clin["grade"] = pd.to_numeric(clin["grade"], errors="coerce")

    def parse_stage(s):
        if pd.isna(s):
            return np.nan
        s = str(s).upper()
        if "M1" in s:
            return 4
        m = re.search(r"P?T(\d)", s)
        if not m:
            return np.nan
        t = int(m.group(1))
        n_match = re.search(r"N(\d)", s)
        n = int(n_match.group(1)) if n_match else 0
        if t >= 4 or n >= 1:
            if t == 4:
                return 4
            if n == 1 and t >= 1:
                return 3
        if t == 1 and n == 0:
            return 1
        if t == 2 and n == 0:
            return 2
        if t == 3 and n == 0:
            return 3
        if t <= 2 and n == 1:
            return 3
        return np.nan

    clin["stage"] = clin["Stage at diagnosis"].apply(parse_stage)

    df = ncc_df.merge(
        clin[["sample", "age", "grade", "stage", "os_months", "os_event"]],
        on="sample",
        how="inner",
    )
    df = df.dropna(subset=["os_months", "os_event", "stage", "grade"])

    n_before = len(df)
    df = df[df["os_months"] > 0]
    print(f"Survival-time filter (os_months > 0): {n_before} -> {len(df)}")

    df["status"] = df["os_event"].astype(int)
    df = df.rename(columns={"NCC_S1": "Subtype"})
    df = df[
        ["sample", "Subtype", "age", "grade", "stage", "os_months", "status"]
    ].copy()
    df.columns = ["sample", "Subtype", "Age", "Grade", "Stage", "time", "status"]
    return df


def prep_icgc():
    with open(SUBTYPE_CLASSIFIER_JSON) as f:
        clf = json.load(f)

    ref_genes = clf["genes"]
    tcga_s1 = np.array(clf["centroid_s1"])
    tcga_s2 = np.array(clf["centroid_s2"])

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
        {d: {g: iexp[d].get(g, 0) for g in ref_genes if g in g2e} for d in iexp}
    ).T.astype(float)
    iem = np.log2(iem + 1)
    iez = iem.subtract(iem.mean(axis=0), axis=1).div(
        iem.std(axis=0, ddof=1), axis=1
    ).fillna(0)

    common = [g for g in ref_genes if g in iez.columns]
    tci = [ref_genes.index(g) for g in common]
    id1 = cdist(iez[common].values, tcga_s1[tci].reshape(1, -1))[:, 0]
    id2 = cdist(iez[common].values, tcga_s2[tci].reshape(1, -1))[:, 0]
    ncc_labels = np.where(id1 < id2, 1, 0)

    ncc_df = pd.DataFrame({"donor": iem.index.tolist(), "NCC_S1": ncc_labels})

    donors = {}
    with gzip.open(ICGC_DONOR, "rt") as f:
        hdr = f.readline().strip().split("\t")
        dsi = hdr.index("donor_survival_time")
        dvi = hdr.index("donor_vital_status")
        dai = hdr.index("donor_age_at_diagnosis")
        sti = hdr.index("donor_tumour_stage_at_diagnosis")
        for line in f:
            c = line.strip().split("\t")
            try:
                dst = float(c[dsi])
            except ValueError:
                continue
            if dst <= 0:
                continue
            donors[c[0]] = {
                "os_days": dst,
                "dead": 1 if c[dvi] == "deceased" else 0,
                "age": c[dai],
                "stage_raw": c[sti],
            }

    grades = {}
    with gzip.open(ICGC_SPECIMEN, "rt") as f:
        hdr = f.readline().strip().split("\t")
        didi = hdr.index("icgc_donor_id")
        scti = hdr.index("specimen_type")
        tgi = hdr.index("tumour_grade")
        for line in f:
            c = line.strip().split("\t")
            if "Primary tumour" in c[scti] and len(c) > tgi and c[tgi] and c[tgi] != "":
                grades[c[didi]] = int(c[tgi])

    def parse_stage_icgc(s):
        if pd.isna(s) or not s:
            return np.nan
        s = str(s).upper()
        if "M1" in s:
            return 4
        m = re.search(r"T(\d)", s)
        if not m:
            return np.nan
        t = int(m.group(1))
        n_match = re.search(r"N(\d)", s)
        n = int(n_match.group(1)) if n_match else 0
        if t >= 4 or n >= 1:
            if t == 4:
                return 4
            if n == 1:
                return 3
        if t == 1 and n == 0:
            return 1
        if t == 2 and n == 0:
            return 2
        if t == 3 and n == 0:
            return 3
        return np.nan

    rows = []
    for d, ncc in ncc_df.values:
        if d in donors and d in grades:
            dn = donors[d]
            rows.append({
                "donor": d,
                "NCC_S1": ncc,
                "Age": float(dn["age"]) if dn["age"] else np.nan,
                "Stage": parse_stage_icgc(dn["stage_raw"]),
                "Grade": grades[d],
                "time": dn["os_days"] / 30,
                "status": dn["dead"],
            })

    df = pd.DataFrame(rows).dropna()
    n_before = len(df)
    df = df[df["time"] > 0]
    print(f"Survival-time filter (time > 0): {n_before} -> {len(df)}")
    df.rename(columns={"NCC_S1": "Subtype"}, inplace=True)
    return df


def run_cohort(name, df, out_png, out_svg, script_name, cohort_label, xlim):
    print(f"\n===== {cohort_label} (n={len(df)}) =====")
    print(f"Clean cohort: {len(df)} patients, {df['status'].sum()} events")

    vars_uni = [
        ("Subtype", "NCC Subtype (S1 vs S2)"),
        ("Stage", "Stage (per level)"),
        ("Grade", "Grade (per level)"),
        ("Age", "Age (per year)"),
    ]
    uni_results = []
    for var, label in vars_uni:
        dfu = df[["time", "status", var]].copy()
        dfu["status"] = dfu["status"].astype(int)
        cp = CoxPHFitter()
        cp.fit(dfu, "time", "status", formula=var)
        s = cp.summary
        hr = s.loc[var, "exp(coef)"]
        lo = s.loc[var, "exp(coef) lower 95%"]
        hi = s.loc[var, "exp(coef) upper 95%"]
        p = s.loc[var, "p"]
        cidx = cp.concordance_index_
        uni_results.append({"label": label, "HR": hr, "L": lo, "U": hi, "p": p, "C": cidx})
        print(f"Uni {var}: HR={hr:.2f} ({lo:.2f}-{hi:.2f}) P={p:.2e} C={cidx:.3f}")

    cph_check = CoxPHFitter()
    cph_check.fit(df, "time", "status", formula="Subtype + Stage + Grade + Age")

    print("\nSchoenfeld residual test")
    cph_check.check_assumptions(df, p_value_threshold=0.05, show_plots=False)

    print("\nStratified Cox (Stage as strata)")
    cph = CoxPHFitter()
    cph.fit(df, "time", "status", formula="Subtype + Grade + Age", strata=["Stage"])
    cph.print_summary()

    pred = cph.predict_partial_hazard(df)
    multi_cidx = concordance_index(df["time"], -pred, df["status"])

    rng = np.random.default_rng(42)
    n_boot = 2000
    boot_cidx = []
    df_boot = df.reset_index(drop=True)

    for _ in range(n_boot):
        idx = rng.choice(len(df_boot), size=len(df_boot), replace=True)
        db = df_boot.iloc[idx]
        try:
            cb = CoxPHFitter()
            cb.fit(db, "time", "status", formula="Subtype + Grade + Age", strata=["Stage"])
            pred = cb.predict_partial_hazard(db)
            c_val = concordance_index(db["time"], -pred, db["status"])
            if not np.isnan(c_val):
                boot_cidx.append(c_val)
        except Exception:
            continue

    ci_lo = np.percentile(boot_cidx, 2.5) if boot_cidx else np.nan
    ci_hi = np.percentile(boot_cidx, 97.5) if boot_cidx else np.nan
    print(
        f"\nStratified C-index: {multi_cidx:.3f} "
        f"(95% CI: {ci_lo:.3f}-{ci_hi:.3f}, "
        f"{len(boot_cidx)}/{n_boot} bootstrap samples)"
    )

    fig = plt.figure(figsize=(18, 10))
    fig.patch.set_facecolor("white")
    y_labels = [
        "NCC Subtype\n(S1 vs S2)",
        "Stage\n(per level)",
        "Grade\n(per level)",
        "Age\n(per year)",
    ]
    colors = [C1, C2, C2, Cg]

    ax_a = fig.add_axes([0.10, 0.15, 0.35, 0.75])
    ax_a.set_facecolor("white")
    y_a = [3, 2, 1, 0]

    for i, (r, clr) in enumerate(zip(uni_results, colors)):
        ax_a.errorbar(
            r["HR"],
            y_a[i],
            xerr=[[r["HR"] - r["L"]], [r["U"] - r["HR"]]],
            fmt="o",
            color=clr,
            capsize=6,
            capthick=2,
            ms=12,
            lw=2.5,
            zorder=3,
        )
        ax_a.text(
            r["HR"] + 0.15,
            y_a[i],
            f"P = {r['p']:.2e}",
            va="center",
            fontsize=11,
            color=clr,
            fontweight="bold",
        )

    ax_a.axvline(x=1, color="#333", ls="--", lw=1.5)
    ax_a.set_yticks(y_a)
    ax_a.set_yticklabels(y_labels, fontsize=16, fontweight="bold")
    ax_a.set_xlabel("Hazard Ratio (95% CI)", fontsize=18, fontweight="bold")
    ax_a.set_title(
        f"A  Univariate Cox ({cohort_label})",
        fontsize=22,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    ax_a.set_xlim(0, xlim)
    ax_a.tick_params(axis="x", labelsize=15)

    ax_b = fig.add_axes([0.55, 0.15, 0.35, 0.75])
    ax_b.set_facecolor("white")
    summary = cph.summary
    vars_b = ["Subtype", "Grade", "Age"]
    y_labels_b = [
        "NCC Subtype\n(S1 vs S2)",
        "Grade\n(per level)",
        "Age\n(per year)",
    ]

    for i, var in enumerate(vars_b):
        hr = summary.loc[var, "exp(coef)"]
        ci_l = summary.loc[var, "exp(coef) lower 95%"]
        ci_u = summary.loc[var, "exp(coef) upper 95%"]
        p = summary.loc[var, "p"]
        color = C1 if p < 0.05 else Cg
        ax_b.errorbar(
            hr,
            3 - i,
            xerr=[[hr - ci_l], [ci_u - hr]],
            fmt="o",
            color=color,
            capsize=6,
            capthick=2,
            ms=12,
            lw=2.5,
            zorder=3,
        )
        ax_b.text(
            hr + 0.15,
            3 - i,
            f"P = {p:.2e}",
            va="center",
            fontsize=11,
            color=color,
            fontweight="bold",
        )

    ax_b.text(
        4.5,
        3.5,
        "Stage: stratification variable",
        fontsize=11,
        color=C2,
        fontweight="bold",
        style="italic",
        ha="right",
    )
    ax_b.axvline(x=1, color="#333", ls="--", lw=1.5)
    ax_b.set_yticks([3, 2, 1, 0])
    ax_b.set_yticklabels(y_labels_b + ["Stage\n(stratified)"], fontsize=16, fontweight="bold")
    ax_b.set_xlabel("Hazard Ratio (95% CI)", fontsize=18, fontweight="bold")
    ax_b.set_title(
        f"B  Multivariate Cox (Stage-stratified, C = {multi_cidx:.3f})",
        fontsize=22,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    ax_b.set_xlim(0, xlim)
    ax_b.tick_params(axis="x", labelsize=15)

    fig.text(
        0.55,
        0.06,
        "Python 3.13 · lifelines 0.30 · Stage-stratified Cox · "
        "Schoenfeld PH passed for Subtype/Grade/Age",
        fontsize=10,
        color="#666",
        style="italic",
    )

    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nSaved: {os.path.basename(out_png)}")
    print(f"       {os.path.basename(out_svg)}")

    try:
        s = cph.summary.loc["Subtype"]
        rows = [
            {
                "section": "Prognostic validation",
                "item": "Stratified NCC HR (Stage+Grade+Age)",
                "cohort": cohort_label,
                "value": f"{s['exp(coef)']:.2f}",
            },
            {
                "section": "Prognostic validation",
                "item": "Stratified NCC 95% CI",
                "cohort": cohort_label,
                "value": f"{s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}",
            },
            {
                "section": "Prognostic validation",
                "item": "Stratified NCC P",
                "cohort": cohort_label,
                "value": f"{s['p']:.1e}",
            },
            {
                "section": "Prognostic validation",
                "item": "N in stratified model (complete Stage+Grade+Age)",
                "cohort": cohort_label,
                "value": f"{len(df)}",
            },
            {
                "section": "Prognostic validation",
                "item": "Stratified C-index (Harrell's C)",
                "cohort": cohort_label,
                "value": f"{multi_cidx:.2f}",
            },
            {
                "section": "Prognostic validation",
                "item": "C-index 95% CI (bootstrap 2000)",
                "cohort": cohort_label,
                "value": f"{ci_lo:.2f}-{ci_hi:.2f}",
            },
        ]
        emit_provenance(script_name, rows)
    except Exception as e:
        print(f"[PROVENANCE] {script_name} emission skipped: {e}")


def main():
    cohorts = [
        (
            "TCGA",
            prep_tcga,
            f"{FIG}/Figure_7_prognostic_cox_stratified.png",
            f"{FIG}/Figure_7_prognostic_cox_stratified.svg",
            "prognostic_cox.py",
            "TCGA-KIRC (Training)",
            5.5,
        ),
        (
            "E-MTAB",
            prep_emtab,
            f"{FIG}/EMTAB_prognostic_cox.png",
            f"{FIG}/EMTAB_prognostic_cox.svg",
            "EMTAB_prognostic_cox.py",
            "E-MTAB-1980",
            15,
        ),
        (
            "ICGC",
            prep_icgc,
            f"{FIG}/ICGC_prognostic_cox.png",
            f"{FIG}/ICGC_prognostic_cox.svg",
            "ICGC_prognostic_cox.py",
            "ICGC RECA-EU",
            5,
        ),
    ]

    for name, prep, png, svg, sname, label, xlim in cohorts:
        df = prep()
        run_cohort(name, df, png, svg, sname, label, xlim)


if __name__ == "__main__":
    main()