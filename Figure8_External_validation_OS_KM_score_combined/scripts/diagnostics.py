"""Diagnostics: proportional-hazards assumption check and survival-time filter audit.

Data preparation is reused from prognostic_cox.py via prep_tcga / prep_emtab /
prep_icgc, so cohort filtering and stage coding stay aligned with the pipeline.

Read-only: fits models, prints to stdout only.
"""

import os
import sys
import warnings

import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import prognostic_cox as pc

TRANSFORMS = ["km", "rank"]


def ph_pvalues(model, df):
    out = {}
    for tt in TRANSFORMS:
        res = proportional_hazard_test(model, df, time_transform=tt)
        for _, row in res.summary.reset_index().iterrows():
            name = row.iloc[0]
            if isinstance(name, tuple):
                name = name[0]
            out.setdefault(str(name), {})[tt] = float(row["p"])
    return out


def report_ph(tag, pvals):
    print(f"  {tag}")
    for cov in sorted(pvals):
        ps = pvals[cov]
        cells = "  ".join(f"{tt}={ps[tt]:.4f}" for tt in TRANSFORMS)
        if all(p < 0.05 for p in ps.values()):
            flag = " <-- VIOLATION"
        elif any(p < 0.05 for p in ps.values()):
            flag = " (transform-dependent)"
        else:
            flag = ""
        print(f"      {cov:<9} {cells}{flag}")


def report_filter(label, t, status=None):
    t = pd.to_numeric(pd.Series(t), errors="coerce")
    n_nan = int(t.isna().sum())
    tv = t.dropna()
    n_neg = int((tv < 0).sum())
    n_zero = int((tv == 0).sum())

    print(f"  {label}")
    print(f"    n = {len(t)}   min = {tv.min():.6g}   max = {tv.max():.6g}")
    print(f"    time < 0 : {n_neg}    time == 0 : {n_zero}    NaN : {n_nan}")

    if status is not None and (n_neg or n_zero):
        s = pd.Series(status).reset_index(drop=True)
        bad = tv.reset_index(drop=True) <= 0
        print(f"    among time<=0, events (status=1): {int(s[bad].sum())}")


def run_ph_check():
    print("=" * 72)
    print("PROPORTIONAL HAZARDS ASSUMPTION CHECK")
    print("Schoenfeld residuals, p < 0.05 in both transforms = violation")
    print("=" * 72)

    cohorts = [
        ("TCGA-KIRC (discovery)", pc.prep_tcga()),
        ("E-MTAB-1980", pc.prep_emtab()),
        ("ICGC RECA-EU", pc.prep_icgc()),
    ]

    for label, df in cohorts:
        print(f"\n{'-' * 72}\n{label}\n{'-' * 72}")
        print(
            f"n = {len(df)}, events = {int(df['status'].sum())}, "
            f"stage levels = {sorted(df['Stage'].unique().tolist())}\n"
        )

        mod_a = CoxPHFitter().fit(
            df, "time", "status", formula="Subtype + Stage + Grade + Age"
        )
        report_ph(
            "[A] DIAGNOSTIC  unstratified, Stage as covariate",
            ph_pvalues(mod_a, df),
        )

        print()
        mod_b = CoxPHFitter().fit(
            df, "time", "status", formula="Subtype + Grade + Age", strata=["Stage"]
        )
        report_ph(
            "[B] REPORTED    stratified by Stage (Stage not testable here)",
            ph_pvalues(mod_b, df),
        )

        s = mod_b.summary
        print(
            f"\n      reported Subtype HR = {s.loc['Subtype', 'exp(coef)']:.2f} "
            f"(P = {s.loc['Subtype', 'p']:.2e})"
        )

    print(f"\n{'=' * 72}\nDONE\n{'=' * 72}")


def run_survtime_audit():
    print("=" * 72)
    print("SURVIVAL-TIME FILTER AUDIT (time > 0 across cohorts)")
    print("=" * 72)

    print("\n[1] TCGA-KIRC  - prognostic_cox.prep_tcga()")
    clin = pd.read_csv(
        os.path.join(pc.CLASSIFIER_DIR, "KIRC_clinical_cleaned.csv"), index_col=0
    )
    sub = pd.read_csv(
        os.path.join(pc.CLASSIFIER_DIR, "subtype_assignment_balanced.csv")
    )
    sub_map = dict(zip(sub["Patient"], sub["Subtype"]))

    clin["Subtype"] = clin.index.map(sub_map)
    clin["Subtype_bin"] = (clin["Subtype"] == "Subtype_1").astype(int)
    clin["Stage_num"] = clin["stage"].map({"I": 1, "II": 2, "III": 3, "IV": 4})
    clin["Grade_num"] = pd.to_numeric(clin["grade"], errors="coerce")
    clin["Age_num"] = pd.to_numeric(clin["age"], errors="coerce")
    clin["os_bin"] = pd.to_numeric(clin["os"], errors="coerce")
    clin["os_t"] = pd.to_numeric(clin["os_time_days"], errors="coerce")

    pre = clin[
        ["os_t", "os_bin", "Subtype_bin", "Stage_num", "Grade_num", "Age_num"]
    ].dropna()
    report_filter("BEFORE time > 0 filter (after dropna)", pre["os_t"], pre["os_bin"])

    post = pc.prep_tcga()
    report_filter("AFTER  filter (model input)", post["time"], post["status"])
    print(f"    -> removed {len(pre) - len(post)} patient(s) by the >0 filter")

    print("\n[2] E-MTAB-1980 - prognostic_cox.prep_emtab()")
    df = pc.prep_emtab()
    report_filter("model input as-is (time in months)", df["time"], df["status"])
    print("    smallest 5 times:", sorted(pd.to_numeric(df["time"]).tolist())[:5])

    print("\n[3] ICGC RECA-EU - prognostic_cox.prep_icgc()")
    df = pc.prep_icgc()
    report_filter("model input as-is (time in months)", df["time"], df["status"])
    print("    smallest 5 times:", sorted(pd.to_numeric(df["time"]).tolist())[:5])

    print("\n[4] TCGA internal consistency")
    t_map = pd.read_csv(pc.SUBTYPE_ASSIGNMENT_CSV)
    t_map = dict(zip(t_map["Patient"], t_map["Subtype"]))

    c = pd.read_csv(pc.TCGA_CLINICAL_CSV, index_col=0)
    c["os_t"] = pd.to_numeric(c["os_time_days"], errors="coerce")
    c["os_b"] = pd.to_numeric(c["os"], errors="coerce")
    c["Subtype"] = c.index.map(t_map)
    c = c.dropna(subset=["os_t", "os_b", "Subtype"])
    c["Subtype_bin"] = (c["Subtype"] == "Subtype_1").astype(int)

    def simple(d, tag):
        s1 = d[d.Subtype_bin == 1]
        s2 = d[d.Subtype_bin == 0]
        lr = logrank_test(s1["os_t"], s2["os_t"], s1["os_b"], s2["os_b"])
        m = CoxPHFitter().fit(
            d[["os_t", "os_b", "Subtype_bin"]], "os_t", "os_b", formula="Subtype_bin"
        )
        st = m.summary.loc["Subtype_bin"]
        print(f"  {tag}")
        print(f"    N={len(d)}  S1={len(s1)}  S2={len(s2)}  events={int(d['os_b'].sum())}")
        print(f"    Log-rank P = {lr.p_value:.3e}")
        print(
            f"    Simple Cox HR = {st['exp(coef)']:.4f} "
            f"({st['exp(coef) lower 95%']:.3f}-{st['exp(coef) upper 95%']:.3f})  "
            f"P = {st['p']:.3e}"
        )

    simple(c, "OLD behaviour (no time > 0 filter) - reference only")
    simple(c[c["os_t"] > 0], "CURRENT (time > 0 filter) - values in summary CSV")
    print("  NOTE: since 2026-08-09 the zero-time patients are excluded everywhere.")
    print("        Univariate (529) and stratified (518) differ only by")
    print("        covariate complete-case (stage / grade / age).")

    print("\n" + "=" * 72)
    print("AUDIT DONE")
    print("=" * 72)


def main():
    run_ph_check()
    run_survtime_audit()
    return 0


if __name__ == "__main__":
    sys.exit(main())