"""Verify all values in NCC_multi_cohort_validation_summary.csv.

Each block mirrors the data processing of its authoritative source script:
    simple_cox_all_cohorts.py     Simple Cox HR, Cox P, log-rank P
    prognostic_cox.py             TCGA stratified Cox, univariate
    EMTAB_prognostic_cox.py       E-MTAB stratified Cox, univariate
    ICGC_prognostic_cox.py        ICGC stratified Cox, univariate
    classifier_concordance_v2.py  68-gene Pearson r
    key14_cross_cohort_v2.py      14-gene direction concordance and Pearson

Cross-cohort effect references use the uniform S1-S2 mean difference on each
cohort's native log scale.
"""

import gzip
import io
import json
import re
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, pearsonr
from statsmodels.stats.multitest import multipletests

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

from config import *

key14 = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36", "SLC27A2",
    "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]

with open(BALANCED_GENES_TXT) as f:
    ncc68 = [line.strip() for line in f]

with open(SUBTYPE_CLASSIFIER_JSON) as f:
    clf = json.load(f)

ref_genes = clf["genes"]
tcga_c1 = np.array(clf["centroid_s1"])
tcga_c2 = np.array(clf["centroid_s2"])


def ok(a, b, tol=0.05):
    return "OK" if abs(a - b) < tol else f"MISMATCH ({a:.4f} vs {b:.4f})"


def ok_str(a, b):
    return "OK" if str(a).strip() == str(b).strip() else f"MISMATCH ({a} vs {b})"


print("=" * 70)
print("FULL TABLE VERIFICATION (unified)")
print("=" * 70)


# --- TCGA-KIRC ---
print("\n--- TCGA-KIRC ---")
tcga = pd.read_csv(TCGA_EXPR, index_col=0)
tcga_sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)

t1 = [c for c in tcga_sub[tcga_sub.Subtype == "Subtype_1"]["Patient"] if c in tcga.columns]
t2 = [c for c in tcga_sub[tcga_sub.Subtype == "Subtype_2"]["Patient"] if c in tcga.columns]

print(f"  S1/S2: {len(t1)}/{len(t2)} (172/361) {ok(len(t1), 172, 1)}")
print(f"  S1%: {len(t1) / (len(t1) + len(t2)) * 100:.1f}% (32.3%)")
print(f"  NCC genes: {len([g for g in ref_genes if g in tcga.index])}/68 (68/68)")

tcga_fc = {
    g: tcga.loc[g, t1].mean() - tcga.loc[g, t2].mean()
    for g in ref_genes
    if g in tcga.index
}
tcga_fc14 = {
    g: tcga.loc[g, t1].mean() - tcga.loc[g, t2].mean()
    for g in key14
    if g in tcga.index
}

clin = pd.read_csv(TCGA_CLINICAL_CSV, index_col=0)
clin["Subtype"] = clin.index.map(dict(zip(tcga_sub["Patient"], tcga_sub["Subtype"])))
clin["os"] = pd.to_numeric(clin["os"], errors="coerce")
clin["os_t"] = pd.to_numeric(clin["os_time_days"], errors="coerce")

cl = clin.dropna(subset=["os", "os_t", "Subtype"])
cl = cl[cl["os_t"] > 0]
print(f"  [Simple Cox] N={len(cl)}, events={int(cl['os'].sum())}")

s1o = cl[cl.Subtype == "Subtype_1"]
s2o = cl[cl.Subtype == "Subtype_2"]
lr = logrank_test(s1o["os_t"], s2o["os_t"], s1o["os"], s2o["os"])
print(f"  Log-rank P: {lr.p_value:.1e} (4.6e-14) {ok(lr.p_value, 4.6e-14, 1e-12)}")

cl["S"] = (cl["Subtype"] == "Subtype_1").astype(int)
cl["status"] = cl["os"].astype(int)
cp1 = CoxPHFitter()
cp1.fit(cl[["os_t", "status", "S"]], "os_t", "status", formula="S")
hr1 = cp1.summary.loc["S", "exp(coef)"]
p1 = cp1.summary.loc["S", "p"]
print(f"  Simple HR: {hr1:.2f} (3.00) {ok(hr1, 3.00, 0.1)}")
print(f"  Simple Cox P: {p1:.1e} (7.0e-13)")

clin["stage_n"] = clin["stage"].map({"I": 1, "II": 2, "III": 3, "IV": 4})
clin["grade_n"] = pd.to_numeric(clin["grade"], errors="coerce")
clin["age_n"] = pd.to_numeric(clin["age"], errors="coerce")

cl2 = clin.dropna(subset=["os", "os_t", "Subtype", "stage_n", "grade_n", "age_n"])
cl2 = cl2[cl2["os_t"] > 0]
cl2["Subtype_bin"] = (cl2["Subtype"] == "Subtype_1").astype(int)
cl2 = cl2[["os_t", "os", "Subtype_bin", "stage_n", "grade_n", "age_n"]].copy()
cl2.columns = ["time", "status", "Subtype", "Stage", "Grade", "Age"]
cl2["status"] = cl2["status"].astype(int)

cp2 = CoxPHFitter()
cp2.fit(cl2, "time", "status", formula="Subtype+Grade+Age", strata=["Stage"])
s = cp2.summary

hr_s = s.loc["Subtype", "exp(coef)"]
lo_s = s.loc["Subtype", "exp(coef) lower 95%"]
hi_s = s.loc["Subtype", "exp(coef) upper 95%"]
p_s = s.loc["Subtype", "p"]
c_s = concordance_index(cl2["time"], -cp2.predict_partial_hazard(cl2), cl2["status"])

print(f"  Strat HR: {hr_s:.2f} ({lo_s:.2f}-{hi_s:.2f}) vs 1.93 (1.39-2.67) {ok(hr_s, 1.93, 0.02)}")
print(f"  Strat P: {p_s:.1e} vs 7.6e-5 {ok(p_s, 7.6e-5, 5e-5)}")
print(f"  Strat C: {c_s:.3f} (Harrell's C) vs 0.720 {ok(c_s, 0.720, 0.02)}")


# --- E-MTAB-1980 ---
print("\n--- E-MTAB-1980 ---")

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
ed1 = cdist(ez[com_em].values, tcga_c1[tci_em].reshape(1, -1))[:, 0]
ed2 = cdist(ez[com_em].values, tcga_c2[tci_em].reshape(1, -1))[:, 0]
e_s1i = np.where(ed1 < ed2)[0]
e_s2i = np.where(ed1 >= ed2)[0]
em_labels = np.where(ed1 < ed2, 1, 0)

print(f"  S1/S2: {len(e_s1i)}/{len(e_s2i)} (26/75) {ok(len(e_s1i), 26, 1)}")
print(f"  NCC genes: {len(com_em)}/68 (66/68) {ok(len(com_em), 66, 1)}")

clin_em_s = pd.read_excel(GEO_SATO_CLINICAL, header=1)
clin_em_s = clin_em_s.rename(
    columns={
        "sample ID": "sample",
        "observation period (month)": "os_m",
        "outcome": "status",
    }
)
clin_em_s["os_m"] = pd.to_numeric(clin_em_s["os_m"], errors="coerce")
clin_em_s["event"] = (clin_em_s["status"].str.lower() == "dead").astype(int)

em_ncc = pd.DataFrame({"sample": list(emat.index), "Subtype": em_labels})
em_df_s = em_ncc.merge(
    clin_em_s[["sample", "os_m", "event"]], on="sample", how="inner"
).dropna(subset=["os_m", "event"])
em_df_s = em_df_s[em_df_s["os_m"] > 0]
em_df_s["status"] = em_df_s["event"].astype(int)

es1 = em_df_s[em_df_s.Subtype == 1]
es2 = em_df_s[em_df_s.Subtype == 0]
lr_e = logrank_test(es1["os_m"], es2["os_m"], es1["status"], es2["status"])

print(f"  [Simple Cox] N={len(em_df_s)}, events={em_df_s['status'].sum()}")
print(
    f"  S1 deaths: {es1['status'].sum()}/{len(es1)} (13/26)  "
    f"S2: {es2['status'].sum()}/{len(es2)} (10/75)  "
    f"{ok_str(es1['status'].sum(), 13)}"
)
print(f"  Log-rank P: {lr_e.p_value:.1e} (3.5e-5) {ok(lr_e.p_value, 3.5e-5, 1e-4)}")

cp_e1 = CoxPHFitter()
cp_e1.fit(em_df_s[["os_m", "status", "Subtype"]], "os_m", "status", formula="Subtype")
se1 = cp_e1.summary.loc["Subtype"]
print(
    f"  Simple HR: {se1['exp(coef)']:.2f} "
    f"({se1['exp(coef) lower 95%']:.2f}-{se1['exp(coef) upper 95%']:.2f}) (4.88) "
    f"{ok(se1['exp(coef)'], 4.88, 0.1)}"
)
print(f"  Simple Cox P: {se1['p']:.1e} (1.8e-4) {ok(se1['p'], 1.8e-4, 1e-4)}")


def parse_stage_emtab(s):
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


clin_em_strat = pd.read_excel(GEO_SATO_CLINICAL, header=1)
clin_em_strat = clin_em_strat.rename(
    columns={
        "sample ID": "sample",
        "Age": "age",
        "Fuhrman grade": "grade",
        "observation period (month)": "os_months",
        "outcome": "status",
    }
)
clin_em_strat["os_months"] = pd.to_numeric(clin_em_strat["os_months"], errors="coerce")
clin_em_strat["os_event"] = (clin_em_strat["status"].str.lower() == "dead").astype(int)
clin_em_strat["age"] = pd.to_numeric(clin_em_strat["age"], errors="coerce")
clin_em_strat["grade"] = pd.to_numeric(clin_em_strat["grade"], errors="coerce")
clin_em_strat["stage"] = clin_em_strat["Stage at diagnosis"].apply(parse_stage_emtab)

em_df_strat = em_ncc.merge(
    clin_em_strat[["sample", "age", "grade", "stage", "os_months", "os_event"]],
    on="sample",
    how="inner",
)
em_df_strat = em_df_strat.dropna(subset=["os_months", "os_event", "stage", "grade"])
em_df_strat["status"] = em_df_strat["os_event"].astype(int)
em_df_strat = em_df_strat.rename(
    columns={"age": "Age", "grade": "Grade", "stage": "Stage", "os_months": "time"}
)
print(f"  [Stratified Cox] N={len(em_df_strat)}, events={em_df_strat['status'].sum()}")

cp_eu = CoxPHFitter()
cp_eu.fit(em_df_strat[["time", "status", "Subtype"]], "time", "status", formula="Subtype")
print(
    f"  Uni HR: {cp_eu.summary.loc['Subtype', 'exp(coef)']:.2f} (4.24) "
    f"{ok(cp_eu.summary.loc['Subtype', 'exp(coef)'], 4.24, 0.1)}"
)

cp_e2 = CoxPHFitter()
cp_e2.fit(em_df_strat, "time", "status", formula="Subtype+Grade+Age", strata=["Stage"])
se2 = cp_e2.summary.loc["Subtype"]
ce2 = concordance_index(
    em_df_strat["time"],
    -cp_e2.predict_partial_hazard(em_df_strat),
    em_df_strat["status"],
)
print(
    f"  Strat HR: {se2['exp(coef)']:.2f} "
    f"({se2['exp(coef) lower 95%']:.2f}-{se2['exp(coef) upper 95%']:.2f}) "
    f"vs 3.37 (1.15-9.83) {ok(se2['exp(coef)'], 3.37, 0.05)}"
)
print(f"  Strat P: {se2['p']:.2e} vs 3.0e-2 {ok(se2['p'], 3.0e-2, 0.01)}")
print(f"  Strat C: {ce2:.3f} (Harrell's C) vs 0.782 {ok(ce2, 0.782, 0.02)}")

em_fc = {g: emat.iloc[e_s1i][g].mean() - emat.iloc[e_s2i][g].mean() for g in com_em}
cm68e = sorted(set(em_fc) & set(tcga_fc))
x68 = [tcga_fc[g] for g in cm68e]
y68 = [em_fc[g] for g in cm68e]
r68e, p68e = pearsonr(x68, y68)
z68e = np.arctanh(r68e)
se68e = 1 / np.sqrt(len(x68) - 3)
print(
    f"  68g r: {r68e:.3f} "
    f"({np.tanh(z68e - 1.96 * se68e):.3f}-{np.tanh(z68e + 1.96 * se68e):.3f}) "
    f"vs 0.881 (0.811-0.925) {ok(r68e, 0.881, 0.02)}"
)
print(f"  68g P: {p68e:.1e} vs 2.0e-22 {ok(p68e, 2.0e-22, 1e-21)}")

em_s1_list = [emat.index[i] for i in range(len(emat)) if ed1[i] < ed2[i]]
em_s2_list = [emat.index[i] for i in range(len(emat)) if ed1[i] >= ed2[i]]

em_14_total = 0
em_14_same = 0
for g in key14:
    if g in emat.columns and g in tcga.index:
        em_14_total += 1
        em_diff = (
            emat.loc[em_s1_list, g].astype(float).mean()
            - emat.loc[em_s2_list, g].astype(float).mean()
        )
        tcga_diff = tcga.loc[g, t1].mean() - tcga.loc[g, t2].mean()
        if (em_diff > 0) == (tcga_diff > 0):
            em_14_same += 1
print(
    f"  14g dir: {em_14_same}/{em_14_total} (10/12) "
    f"{ok_str(em_14_same, 10)}  [total={em_14_total}]"
)

em_fc14 = {
    g: emat.iloc[e_s1i][g].mean() - emat.iloc[e_s2i][g].mean()
    for g in key14
    if g in emat.columns
}
cm14e = sorted(set(tcga_fc14) & set(em_fc14))
x14 = [tcga_fc14[g] for g in cm14e]
y14 = [em_fc14[g] for g in cm14e]
r14e, p14e = pearsonr(x14, y14)
z14e = np.arctanh(r14e)
se14e = 1 / np.sqrt(len(x14) - 3)
print(
    f"  14g r: {r14e:.3f} "
    f"({np.tanh(z14e - 1.96 * se14e):.3f}-{np.tanh(z14e + 1.96 * se14e):.3f}) "
    f"vs 0.940 (0.795-0.983) {ok(r14e, 0.940, 0.02)}"
)
print(f"  14g P: {p14e:.1e} vs 5.5e-6 {ok(p14e, 5.5e-6, 1e-5)}")


# --- ICGC RECA-EU ---
print("\n--- ICGC RECA-EU ---")

all_sym = list(set(ref_genes + key14))
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
id1 = cdist(iez[com_ic].values, tcga_c1[tci_ic].reshape(1, -1))[:, 0]
id2 = cdist(iez[com_ic].values, tcga_c2[tci_ic].reshape(1, -1))[:, 0]
i_s1i = np.where(id1 < id2)[0]
i_s2i = np.where(id1 >= id2)[0]
ic_labels = np.where(id1 < id2, 1, 0)

print(f"  S1/S2: {len(i_s1i)}/{len(i_s2i)} (21/70) {ok(len(i_s1i), 21, 1)}")
print(f"  NCC genes: {len(com_ic)}/68 (68/68) {ok(len(com_ic), 68, 1)}")

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
            ic_os[c[0]] = {
                "time": ds / 30,
                "dead": 1 if c[dvi] == "deceased" else 0,
            }

ic_rows = []
for i, d in enumerate(iem.index):
    if d in ic_os:
        ic_rows.append(
            {
                "donor": d,
                "Subtype": ic_labels[i],
                "time": ic_os[d]["time"],
                "status": ic_os[d]["dead"],
            }
        )

ic_df_s = pd.DataFrame(ic_rows).dropna(subset=["time", "status"])
ic_df_s["status"] = ic_df_s["status"].astype(int)
ic_df_s = ic_df_s[ic_df_s["time"] > 0]

is1 = ic_df_s[ic_df_s.Subtype == 1]
is2 = ic_df_s[ic_df_s.Subtype == 0]
lr_i = logrank_test(is1["time"], is2["time"], is1["status"], is2["status"])

print(f"  [Simple Cox] N={len(ic_df_s)}, events={ic_df_s['status'].sum()}")
print(f"  Log-rank P: {lr_i.p_value:.2f} (0.30) {ok(lr_i.p_value, 0.30, 0.1)}")

cp_i1 = CoxPHFitter()
cp_i1.fit(ic_df_s[["time", "status", "Subtype"]], "time", "status", formula="Subtype")
si1 = cp_i1.summary.loc["Subtype"]
print(
    f"  Simple HR: {si1['exp(coef)']:.2f} "
    f"({si1['exp(coef) lower 95%']:.2f}-{si1['exp(coef) upper 95%']:.2f}) (1.51) "
    f"{ok(si1['exp(coef)'], 1.51, 0.1)}"
)
print(f"  Simple Cox P: {si1['p']:.2f} (0.30) {ok(si1['p'], 0.30, 0.1)}")

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

ic_grade = {}
with gzip.open(ICGC_SPECIMEN, "rt") as f:
    hdr = f.readline().strip().split("\t")
    didi = hdr.index("icgc_donor_id")
    scti = hdr.index("specimen_type")
    tgi = hdr.index("tumour_grade")
    for line in f:
        c = line.strip().split("\t")
        if "Primary tumour" in c[scti] and len(c) > tgi and c[tgi] and c[tgi] != "":
            ic_grade[c[didi]] = int(c[tgi])


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


ic_ncc_df = pd.DataFrame({"donor": iem.index.tolist(), "NCC_S1": ic_labels})
ic_strat_rows = []
for d, ncc in ic_ncc_df.values:
    if d in donors and d in ic_grade:
        dn = donors[d]
        ic_strat_rows.append({
            "donor": d,
            "NCC_S1": ncc,
            "Age": float(dn["age"]) if dn["age"] else np.nan,
            "Stage": parse_stage_icgc(dn["stage_raw"]),
            "Grade": ic_grade[d],
            "time": dn["os_days"] / 30,
            "status": dn["dead"],
        })

ic_df_strat = pd.DataFrame(ic_strat_rows).dropna()
ic_df_strat = ic_df_strat.rename(columns={"NCC_S1": "Subtype"})
print(f"  [Stratified Cox] N={len(ic_df_strat)}, events={ic_df_strat['status'].sum()}")

cp_iu = CoxPHFitter()
cp_iu.fit(
    ic_df_strat[["time", "status", "Subtype"]],
    "time",
    "status",
    formula="Subtype",
)
print(
    f"  Uni HR: {cp_iu.summary.loc['Subtype', 'exp(coef)']:.2f} (1.35) "
    f"{ok(cp_iu.summary.loc['Subtype', 'exp(coef)'], 1.35, 0.1)}"
)

cp_i2 = CoxPHFitter()
cp_i2.fit(ic_df_strat, "time", "status", formula="Subtype+Grade+Age", strata=["Stage"])
si2 = cp_i2.summary.loc["Subtype"]
ci2 = concordance_index(
    ic_df_strat["time"],
    -cp_i2.predict_partial_hazard(ic_df_strat),
    ic_df_strat["status"],
)
print(
    f"  Strat HR: {si2['exp(coef)']:.2f} "
    f"({si2['exp(coef) lower 95%']:.2f}-{si2['exp(coef) upper 95%']:.2f}) "
    f"vs 1.28 (0.55-2.98) {ok(si2['exp(coef)'], 1.28, 0.05)}"
)
print(f"  Strat P: {si2['p']:.2f} vs 0.57 {ok(si2['p'], 0.57, 0.05)}")
print(f"  Strat C: {ci2:.3f} (Harrell's C) vs 0.673 {ok(ci2, 0.673, 0.02)}")

ic_fc = {g: iem.iloc[i_s1i][g].mean() - iem.iloc[i_s2i][g].mean() for g in com_ic}
cm68i = sorted(set(ic_fc) & set(tcga_fc))
x68i = [tcga_fc[g] for g in cm68i]
y68i = [ic_fc[g] for g in cm68i]
r68i, p68i = pearsonr(x68i, y68i)
z68i = np.arctanh(r68i)
se68i = 1 / np.sqrt(len(x68i) - 3)
print(
    f"  68g r: {r68i:.3f} "
    f"({np.tanh(z68i - 1.96 * se68i):.3f}-{np.tanh(z68i + 1.96 * se68i):.3f}) "
    f"vs 0.735 (0.602-0.828) {ok(r68i, 0.735, 0.02)}"
)
print(f"  68g P: {p68i:.1e} vs 1.0e-12 {ok(p68i, 1.0e-12, 1e-11)}")

ic_s1_list = [iem.index[i] for i in range(len(iem)) if id1[i] < id2[i]]
ic_s2_list = [iem.index[i] for i in range(len(iem)) if id1[i] >= id2[i]]

ic_14_total = 0
ic_14_same = 0
for g in key14:
    if g in iem.columns and g in tcga.index:
        ic_14_total += 1
        ic_diff = (
            iem.loc[ic_s1_list, g].astype(float).mean()
            - iem.loc[ic_s2_list, g].astype(float).mean()
        )
        tcga_diff = tcga.loc[g, t1].mean() - tcga.loc[g, t2].mean()
        if (ic_diff > 0) == (tcga_diff > 0):
            ic_14_same += 1
print(
    f"  14g dir: {ic_14_same}/{ic_14_total} (12/14) "
    f"{ok_str(ic_14_same, 12)}  [total={ic_14_total}]"
)

ic_fc14 = {
    g: iem.iloc[i_s1i][g].mean() - iem.iloc[i_s2i][g].mean()
    for g in key14
    if g in iem.columns
}
cm14i = sorted(set(tcga_fc14) & set(ic_fc14))
x14i = [tcga_fc14[g] for g in cm14i]
y14i = [ic_fc14[g] for g in cm14i]
r14i, p14i = pearsonr(x14i, y14i)
z14i = np.arctanh(r14i)
se14i = 1 / np.sqrt(len(x14i) - 3)
print(
    f"  14g r: {r14i:.3f} "
    f"({np.tanh(z14i - 1.96 * se14i):.3f}-{np.tanh(z14i + 1.96 * se14i):.3f}) "
    f"vs 0.740 (0.344-0.912) {ok(r14i, 0.740, 0.02)}"
)
print(f"  14g P: {p14i:.1e} vs 2.5e-3 {ok(p14i, 2.5e-3, 1e-3)}")


# --- CPTAC ccRCC ---
print("\n--- CPTAC ccRCC ---")

rnaseq_cptac = pd.read_csv(CPTAC_RNA, sep="\t", index_col=0)
rnaz_cp = rnaseq_cptac.subtract(rnaseq_cptac.mean(axis=1), axis=0).div(
    rnaseq_cptac.std(axis=1, ddof=1), axis=0
).fillna(0)

avail_clf = [g for g in ref_genes if g in rnaz_cp.index]
tci_cpv = [ref_genes.index(g) for g in avail_clf]
expr_clf_z = rnaz_cp.loc[avail_clf].values.T
cc1 = np.array(clf["centroid_s1"])[tci_cpv]
cc2 = np.array(clf["centroid_s2"])[tci_cpv]
cd1 = cdist(expr_clf_z, cc1.reshape(1, -1))[:, 0]
cd2 = cdist(expr_clf_z, cc2.reshape(1, -1))[:, 0]
cp_s1s = np.array(rnaseq_cptac.columns)[cd1 < cd2]
cp_s2s = np.array(rnaseq_cptac.columns)[cd1 >= cd2]

print(f"  NCC S1/S2: {len(cp_s1s)}/{len(cp_s2s)}")
print(f"  NCC genes: {len(avail_clf)}/68")

cp_fc68 = {
    g: rnaseq_cptac.loc[g, cd1 < cd2].mean() - rnaseq_cptac.loc[g, cd1 >= cd2].mean()
    for g in avail_clf
}
cm68cp = sorted(set(tcga_fc) & set(cp_fc68))
x68cp = [tcga_fc[g] for g in cm68cp]
y68cp = [cp_fc68[g] for g in cm68cp]
r68cp, p68cp = pearsonr(x68cp, y68cp)
z68cp = np.arctanh(r68cp)
se68cp = 1 / np.sqrt(len(x68cp) - 3)
print(
    f"  68g r: {r68cp:.3f} "
    f"({np.tanh(z68cp - 1.96 * se68cp):.3f}-{np.tanh(z68cp + 1.96 * se68cp):.3f}) "
    f"vs 0.928 (0.885-0.956) {ok(r68cp, 0.928, 0.02)}"
)
print(f"  68g P: {p68cp:.1e} vs 8.9e-29 {ok(p68cp, 8.9e-29, 1e-10)}")

prot_tumor = pd.read_csv(CPTAC_PROTEIN, sep="\t", index_col=0)
f_avail = [g for g in key14 if g in prot_tumor.index]

f_data = []
for g in f_avail:
    s1 = prot_tumor.loc[
        g, [c for c in cp_s1s if c in prot_tumor.columns]
    ].values.astype(float)
    s2 = prot_tumor.loc[
        g, [c for c in cp_s2s if c in prot_tumor.columns]
    ].values.astype(float)
    s1 = s1[~np.isnan(s1)]
    s2 = s2[~np.isnan(s2)]
    if len(s1) < 3 or len(s2) < 3:
        continue
    fc = s1.mean() - s2.mean()
    _, p = mannwhitneyu(s1, s2)
    f_data.append({"gene": g, "log2FC": fc, "p": p})

f_df = pd.DataFrame(f_data)
if len(f_df) > 0:
    _, f_df["p_adj"], _, _ = multipletests(f_df["p"].values, method="fdr_bh")
else:
    f_df["p_adj"] = f_df["p"]

valid_genes = [g for g in f_df["gene"] if g in tcga_fc14]
agreement = sum(
    (f_df.loc[f_df.gene == g, "log2FC"].values[0] > 0) == (tcga_fc14[g] > 0)
    for g in valid_genes
)
n_bh = int((f_df["p_adj"] < 0.05).sum())
print(
    f"  14g dir: {agreement}/{len(valid_genes)} (9/12) "
    f"{ok_str(agreement, 9)}  [total={len(f_df)}]"
)
print(f"  14g BH-sig: {n_bh}/{len(f_df)} (5/12) {ok_str(n_bh, 5)}")

cp_prot_fc = {}
for g in key14:
    if g in prot_tumor.index:
        v1 = prot_tumor.loc[
            g, [c for c in cp_s1s if c in prot_tumor.columns]
        ].astype(float)
        v2 = prot_tumor.loc[
            g, [c for c in cp_s2s if c in prot_tumor.columns]
        ].astype(float)
        v1 = v1[~np.isnan(v1)]
        v2 = v2[~np.isnan(v2)]
        if len(v1) >= 3 and len(v2) >= 3:
            cp_prot_fc[g] = v1.mean() - v2.mean()

cm14cp = sorted(set(tcga_fc14) & set(cp_prot_fc))
x14cp = [tcga_fc14[g] for g in cm14cp]
y14cp = [cp_prot_fc[g] for g in cm14cp]
r14cp, p14cp = pearsonr(x14cp, y14cp)
z14cp = np.arctanh(r14cp)
se14cp = 1 / np.sqrt(len(x14cp) - 3)
print(
    f"  14g r: {r14cp:.3f} "
    f"({np.tanh(z14cp - 1.96 * se14cp):.3f}-{np.tanh(z14cp + 1.96 * se14cp):.3f}) "
    f"vs 0.848 (0.535-0.957) {ok(r14cp, 0.848, 0.02)}"
)
print(f"  14g P: {p14cp:.1e} vs 4.9e-4 {ok(p14cp, 4.9e-4, 1e-4)}")


print("\n" + "=" * 70)
print("VERIFICATION COMPLETE")
print("Any MISMATCH above indicates a value that needs correction")
print("=" * 70)