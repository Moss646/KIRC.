"""Simple Cox + log-rank across three cohorts with OS data.

Generates NCC HR, Cox P, log-rank P for TCGA-KIRC, E-MTAB-1980, ICGC RECA-EU.
"""

import gzip
import io
import json
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy.spatial.distance import cdist

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

from config import *

with open(SUBTYPE_CLASSIFIER_JSON) as f:
    clf = json.load(f)

ref_genes = clf["genes"]
c1 = np.array(clf["centroid_s1"])
c2 = np.array(clf["centroid_s2"])

print("=" * 60)
print("SIMPLE COX + LOG-RANK - ALL COHORTS")
print("=" * 60)


# --- TCGA-KIRC ---
print("\n[TCGA-KIRC]")
tcga_expr = pd.read_csv(TCGA_EXPR, index_col=0)
tcga_sub = pd.read_csv(SUBTYPE_ASSIGNMENT_CSV)
t_map = dict(zip(tcga_sub["Patient"], tcga_sub["Subtype"]))

t_clin = pd.read_csv(TCGA_CLINICAL_CSV, index_col=0)
t_clin["os_t"] = pd.to_numeric(t_clin["os_time_days"], errors="coerce")
t_clin["os_b"] = pd.to_numeric(t_clin["os"], errors="coerce")
t_clin["Subtype"] = t_clin.index.map(t_map)
t_clin = t_clin.dropna(subset=["os_t", "os_b", "Subtype"])

n_before = len(t_clin)
t_clin = t_clin[t_clin["os_t"] > 0]
print(f"  Survival-time filter (os_t > 0): {n_before} -> {len(t_clin)}")

t_clin["Subtype_bin"] = (t_clin["Subtype"] == "Subtype_1").astype(int)
t_s1 = t_clin[t_clin.Subtype_bin == 1]
t_s2 = t_clin[t_clin.Subtype_bin == 0]

lr_t = logrank_test(t_s1["os_t"], t_s2["os_t"], t_s1["os_b"], t_s2["os_b"])

cph_t = CoxPHFitter()
cph_t.fit(t_clin[["os_t", "os_b", "Subtype_bin"]], "os_t", "os_b", formula="Subtype_bin")
st = cph_t.summary.loc["Subtype_bin"]

print(f"  S1={len(t_s1)}, S2={len(t_s2)}, events={int(t_clin['os_b'].sum())}")
print(f"  Log-rank P = {lr_t.p_value:.1e}")
print(
    f"  Simple Cox HR = {st['exp(coef)']:.2f} "
    f"(95% CI: {st['exp(coef) lower 95%']:.2f}-{st['exp(coef) upper 95%']:.2f})  "
    f"P = {st['p']:.1e}"
)


# --- E-MTAB-1980 ---
print("\n[E-MTAB-1980]")
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
em_df = em_ncc.merge(
    clin_em[["sample", "os_m", "event"]], on="sample", how="inner"
).dropna(subset=["os_m", "event"])

n_before = len(em_df)
em_df = em_df[em_df["os_m"] > 0]
print(f"  Survival-time filter (os_m > 0): {n_before} -> {len(em_df)}")

em_df["status"] = em_df["event"].astype(int)
em_s1 = em_df[em_df.Subtype == 1]
em_s2 = em_df[em_df.Subtype == 0]

lr_e = logrank_test(em_s1["os_m"], em_s2["os_m"], em_s1["status"], em_s2["status"])

cph_e = CoxPHFitter()
cph_e.fit(em_df[["os_m", "status", "Subtype"]], "os_m", "status", formula="Subtype")
se = cph_e.summary.loc["Subtype"]

print(f"  S1={len(em_s1)}, S2={len(em_s2)}, events={int(em_df['status'].sum())}")
print(f"  Log-rank P = {lr_e.p_value:.1e}")
print(
    f"  Simple Cox HR = {se['exp(coef)']:.2f} "
    f"(95% CI: {se['exp(coef) lower 95%']:.2f}-{se['exp(coef) upper 95%']:.2f})  "
    f"P = {se['p']:.1e}"
)


# --- ICGC RECA-EU ---
print("\n[ICGC RECA-EU]")
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
            {
                "donor": d,
                "Subtype": ic_labels[i],
                "time": ic_os[d]["time"],
                "status": ic_os[d]["dead"],
            }
        )

ic_df = pd.DataFrame(ic_rows).dropna(subset=["time", "status"])
ic_df["status"] = ic_df["status"].astype(int)

n_before = len(ic_df)
ic_df = ic_df[ic_df["time"] > 0]
print(f"  Survival-time filter (time > 0): {n_before} -> {len(ic_df)}")

ic_s1 = ic_df[ic_df.Subtype == 1]
ic_s2 = ic_df[ic_df.Subtype == 0]

lr_i = logrank_test(ic_s1["time"], ic_s2["time"], ic_s1["status"], ic_s2["status"])

cph_i = CoxPHFitter()
cph_i.fit(ic_df[["time", "status", "Subtype"]], "time", "status", formula="Subtype")
si = cph_i.summary.loc["Subtype"]

print(f"  S1={len(ic_s1)}, S2={len(ic_s2)}, events={int(ic_df['status'].sum())}")
print(f"  Log-rank P = {lr_i.p_value:.2e}")
print(
    f"  Simple Cox HR = {si['exp(coef)']:.2f} "
    f"(95% CI: {si['exp(coef) lower 95%']:.2f}-{si['exp(coef) upper 95%']:.2f})  "
    f"P = {si['p']:.2e}"
)


# --- Summary ---
thr = st["exp(coef)"]
tp = st["p"]
tlp = lr_t.p_value
shr = se["exp(coef)"]
sp = se["p"]
slp = lr_e.p_value
ihr = si["exp(coef)"]
ipv = si["p"]
ilp = lr_i.p_value

print("\n" + "=" * 60)
print("SUMMARY - TABLE VALUES")
print("=" * 60)
print(f"TCGA:     NCC HR={thr:.2f}   Cox P={tp:.1e}   Log-rank P={tlp:.1e}")
print(f"E-MTAB:   NCC HR={shr:.2f}   Cox P={sp:.1e}   Log-rank P={slp:.1e}")
print(f"ICGC:     NCC HR={ihr:.2f}   Cox P={ipv:.2e}   Log-rank P={ilp:.2e}")

print("\nConsistent with table?")
ok_t = "OK" if abs(thr - 3.00) < 0.1 else "MISMATCH"
ok_e = "OK" if abs(shr - 4.88) < 0.1 else "MISMATCH"
ok_i = "OK" if abs(ihr - 1.51) < 0.1 else "MISMATCH"
print(f"  TCGA:    {ok_t} (computed={thr:.2f}, table=3.00)")
print(f"  E-MTAB:  {ok_e} (computed={shr:.2f}, table=4.88)")
print(f"  ICGC:    {ok_i} (computed={ihr:.2f}, table=1.51)")

try:
    rows = [
        {"section": "Prognostic validation", "item": "NCC HR",
         "cohort": "TCGA-KIRC (Training)", "value": f"{thr:.2f}"},
        {"section": "Prognostic validation", "item": "Cox P",
         "cohort": "TCGA-KIRC (Training)", "value": f"{tp:.1e}"},
        {"section": "Prognostic validation", "item": "Log-rank P",
         "cohort": "TCGA-KIRC (Training)", "value": f"{tlp:.1e}"},

        {"section": "Prognostic validation", "item": "NCC HR",
         "cohort": "E-MTAB-1980", "value": f"{shr:.2f}"},
        {"section": "Prognostic validation", "item": "Cox P",
         "cohort": "E-MTAB-1980", "value": f"{sp:.1e}"},
        {"section": "Prognostic validation", "item": "Log-rank P",
         "cohort": "E-MTAB-1980", "value": f"{slp:.1e}"},

        {"section": "Prognostic validation", "item": "NCC HR",
         "cohort": "ICGC RECA-EU", "value": f"{ihr:.2f}"},
        {"section": "Prognostic validation", "item": "Cox P",
         "cohort": "ICGC RECA-EU", "value": f"{ipv:.1e}"},
        {"section": "Prognostic validation", "item": "Log-rank P",
         "cohort": "ICGC RECA-EU", "value": f"{ilp:.1e}"},
    ]
    emit_provenance("simple_cox_all_cohorts.py", rows)
except Exception as e:
    print(f"[PROVENANCE] simple_cox_all_cohorts emission skipped: {e}")