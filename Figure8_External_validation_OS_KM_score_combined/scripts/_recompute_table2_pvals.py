"""Recompute Table 2 P values.

Log-rank, univariable Cox, stage-stratified Cox, and C-index for the three
cohorts. Mirrors simple_cox_all_cohorts.py and prognostic_cox.py exactly
(seed 42, n_boot 2000). Writes results to _recompute_table2_pvals.txt.
"""

import gzip
import json
import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from scipy.spatial.distance import cdist

from config import *
from prognostic_cox import prep_tcga, prep_emtab, prep_icgc

warnings.filterwarnings("ignore")

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "_recompute_table2_pvals.txt"
)

_lines = []


def say(s=""):
    _lines.append(str(s))
    print(s)


with open(SUBTYPE_CLASSIFIER_JSON) as f:
    clf = json.load(f)
ref_genes = clf["genes"]
c1 = np.array(clf["centroid_s1"])
c2 = np.array(clf["centroid_s2"])


def fmt_p(p):
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


say("=" * 70)
say("PART A  Log-rank P + univariable Cox (all patients with OS data)")
say("=" * 70)

say("\n[TCGA-KIRC]")
t_clin = pd.read_csv(TCGA_CLINICAL_CSV, index_col=0)
t_clin["os_t"] = pd.to_numeric(t_clin["os_time_days"], errors="coerce")
t_clin["os_b"] = pd.to_numeric(t_clin["os"], errors="coerce")

tcga_sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)
t_map = dict(zip(tcga_sub["Patient"], tcga_sub["Subtype"]))
t_clin["Subtype"] = t_clin.index.map(t_map)
t_clin = t_clin.dropna(subset=["os_t", "os_b", "Subtype"])

n0 = len(t_clin)
t_clin = t_clin[t_clin["os_t"] > 0]
say(f"  Survival-time filter (os_t > 0): {n0} -> {len(t_clin)}")

t_clin["Subtype_bin"] = (t_clin["Subtype"] == "Subtype_1").astype(int)
t_s1 = t_clin[t_clin.Subtype_bin == 1]
t_s2 = t_clin[t_clin.Subtype_bin == 0]

lr_t = logrank_test(t_s1["os_t"], t_s2["os_t"], t_s1["os_b"], t_s2["os_b"])

cph_t = CoxPHFitter()
cph_t.fit(t_clin[["os_t", "os_b", "Subtype_bin"]], "os_t", "os_b", formula="Subtype_bin")
st = cph_t.summary.loc["Subtype_bin"]

say(f"  S1={len(t_s1)}, S2={len(t_s2)}, events={int(t_clin['os_b'].sum())}")
say(f"  Log-rank P = {lr_t.p_value:.6g}  -> table fmt: {fmt_p(lr_t.p_value)}")
say(
    f"  Univariable HR = {st['exp(coef)']:.3f} "
    f"(95% CI {st['exp(coef) lower 95%']:.3f}-{st['exp(coef) upper 95%']:.3f})  "
    f"P = {fmt_p(st['p'])}"
)

say("\n[E-MTAB-1980]")
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
ez = emat.subtract(emat.mean(axis=0), axis=1).div(emat.std(axis=0, ddof=1), axis=1).fillna(0)

com_em = [g for g in ref_genes if g in ez.columns]
tci_em = [ref_genes.index(g) for g in com_em]
ed1 = cdist(ez[com_em].values, c1[tci_em].reshape(1, -1))[:, 0]
ed2 = cdist(ez[com_em].values, c2[tci_em].reshape(1, -1))[:, 0]
em_labels = np.where(ed1 < ed2, 1, 0)

clin_em = pd.read_excel(GEO_SATO_CLINICAL, header=1)
clin_em = clin_em.rename(
    columns={
        "sample ID": "sample",
        "observation period (month)": "os_m",
        "outcome": "status",
    }
)
clin_em["os_m"] = pd.to_numeric(clin_em["os_m"], errors="coerce")
clin_em["event"] = (clin_em["status"].str.lower() == "dead").astype(int)

em_ncc = pd.DataFrame({"sample": list(emat.index), "Subtype": em_labels})
em_df = em_ncc.merge(clin_em[["sample", "os_m", "event"]], on="sample", how="inner")
em_df = em_df.dropna(subset=["os_m", "event"])

n0 = len(em_df)
em_df = em_df[em_df["os_m"] > 0]
say(f"  Survival-time filter (os_m > 0): {n0} -> {len(em_df)}")

em_df["status"] = em_df["event"].astype(int)
em_s1 = em_df[em_df.Subtype == 1]
em_s2 = em_df[em_df.Subtype == 0]

lr_e = logrank_test(em_s1["os_m"], em_s2["os_m"], em_s1["status"], em_s2["status"])

cph_e = CoxPHFitter()
cph_e.fit(em_df[["os_m", "status", "Subtype"]], "os_m", "status", formula="Subtype")
se = cph_e.summary.loc["Subtype"]

say(f"  S1-like={len(em_s1)}, S2-like={len(em_s2)}, events={int(em_df['status'].sum())}")
say(f"  Log-rank P = {lr_e.p_value:.6g}  -> table fmt: {fmt_p(lr_e.p_value)}")
say(
    f"  Univariable HR = {se['exp(coef)']:.3f} "
    f"(95% CI {se['exp(coef) lower 95%']:.3f}-{se['exp(coef) upper 95%']:.3f})  "
    f"P = {fmt_p(se['p'])}"
)

say("\n[ICGC RECA-EU]")
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
iez = iem.subtract(iem.mean(axis=0), axis=1).div(iem.std(axis=0, ddof=1), axis=1).fillna(0)

com_ic = [g for g in ref_genes if g in iez.columns]
tci_ic = [ref_genes.index(g) for g in com_ic]
id1 = cdist(iez[com_ic].values, c1[tci_ic].reshape(1, -1))[:, 0]
id2 = cdist(iez[com_ic].values, c2[tci_ic].reshape(1, -1))[:, 0]
ic_labels = np.where(id1 < id2, 1, 0)

ic_os = {}
with gzip.open(ICGC_DONOR, "rt") as f:
    hdr = f.readline().strip().split("\t")
    dsi = hdr.index("donor_survival_time")
    dvi = hdr.index("donor_vital_status")
    for line in f:
        c = line.strip().split("\t")
        try:
            ds = float(c[dsi])
        except ValueError:
            continue
        if ds > 0:
            ic_os[c[0]] = {"time": ds / 30, "dead": 1 if c[dvi] == "deceased" else 0}

ic_rows = []
for i, d in enumerate(iem.index):
    if d in ic_os:
        ic_rows.append(
            {"donor": d, "Subtype": ic_labels[i], "time": ic_os[d]["time"], "status": ic_os[d]["dead"]}
        )

ic_df = pd.DataFrame(ic_rows).dropna(subset=["time", "status"])
ic_df["status"] = ic_df["status"].astype(int)

n0 = len(ic_df)
ic_df = ic_df[ic_df["time"] > 0]
say(f"  Survival-time filter (time > 0): {n0} -> {len(ic_df)}")

ic_s1 = ic_df[ic_df.Subtype == 1]
ic_s2 = ic_df[ic_df.Subtype == 0]

lr_i = logrank_test(ic_s1["time"], ic_s2["time"], ic_s1["status"], ic_s2["status"])

cph_i = CoxPHFitter()
cph_i.fit(ic_df[["time", "status", "Subtype"]], "time", "status", formula="Subtype")
si = cph_i.summary.loc["Subtype"]

say(f"  S1-like={len(ic_s1)}, S2-like={len(ic_s2)}, events={int(ic_df['status'].sum())}")
say(f"  Log-rank P = {lr_i.p_value:.6g}  -> table fmt: {fmt_p(lr_i.p_value)}")
say(
    f"  Univariable HR = {si['exp(coef)']:.3f} "
    f"(95% CI {si['exp(coef) lower 95%']:.3f}-{si['exp(coef) upper 95%']:.3f})  "
    f"P = {fmt_p(si['p'])}"
)

say("\n" + "=" * 70)
say("PART B  Stage-stratified Cox + C-index (complete Stage+Grade+Age cases)")
say("=" * 70)

cohorts = [
    ("TCGA-KIRC", prep_tcga()),
    ("E-MTAB-1980", prep_emtab()),
    ("ICGC RECA-EU", prep_icgc()),
]

for label, df in cohorts:
    say(f"\n[{label}] n={len(df)}, events={int(df['status'].sum())}")

    cph = CoxPHFitter()
    cph.fit(df, "time", "status", formula="Subtype + Grade + Age", strata=["Stage"])
    s = cph.summary.loc["Subtype"]

    pred = cph.predict_partial_hazard(df)
    multi_cidx = concordance_index(df["time"], -pred, df["status"])

    say(
        f"  Stage-stratified HR = {s['exp(coef)']:.3f} "
        f"(95% CI {s['exp(coef) lower 95%']:.3f}-{s['exp(coef) upper 95%']:.3f})  "
        f"P = {s['p']:.6g} -> table fmt: {fmt_p(s['p'])}"
    )

    rng = np.random.default_rng(42)
    n_boot = 2000
    df_boot = df.reset_index(drop=True)
    boot_cidx = []

    for _ in range(n_boot):
        idx = rng.choice(len(df_boot), size=len(df_boot), replace=True)
        db = df_boot.iloc[idx]
        try:
            cb = CoxPHFitter()
            cb.fit(db, "time", "status", formula="Subtype + Grade + Age", strata=["Stage"])
            pred_b = cb.predict_partial_hazard(db)
            c_val = concordance_index(db["time"], -pred_b, db["status"])
            if not np.isnan(c_val):
                boot_cidx.append(c_val)
        except Exception:
            continue

    ci_lo = np.percentile(boot_cidx, 2.5)
    ci_hi = np.percentile(boot_cidx, 97.5)
    say(
        f"  C-index = {multi_cidx:.3f} "
        f"(95% CI {ci_lo:.3f}-{ci_hi:.3f}, {len(boot_cidx)}/{n_boot} bootstrap)"
    )

say("\n" + "=" * 70)
say("SUMMARY vs Table 2")
say("=" * 70)
say(
    f"TCGA-KIRC : Log-rank P = {fmt_p(lr_t.p_value)} | "
    f"Uni HR {st['exp(coef)']:.2f} "
    f"({st['exp(coef) lower 95%']:.2f}-{st['exp(coef) upper 95%']:.2f})"
)
say(
    f"E-MTAB    : Log-rank P = {fmt_p(lr_e.p_value)} | "
    f"Uni HR {se['exp(coef)']:.2f} "
    f"({se['exp(coef) lower 95%']:.2f}-{se['exp(coef) upper 95%']:.2f})"
)
say(
    f"ICGC      : Log-rank P = {fmt_p(lr_i.p_value)} | "
    f"Uni HR {si['exp(coef)']:.2f} "
    f"({si['exp(coef) lower 95%']:.2f}-{si['exp(coef) upper 95%']:.2f})"
)
say("(Stage-stratified HR / P / C-index values printed in PART B above)")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(_lines))

print(f"\nSaved -> {OUT}")