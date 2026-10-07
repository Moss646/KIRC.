"""KM survival curves for external validation cohorts, publication size.

Figure is 180 x 90 mm at 600 dpi. Analysis is identical to
external_validation_KM.py; only the figure block differs.
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
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.spatial.distance import cdist

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.linewidth"] = 0.8

from config import *

FIG = FIGURES_DIR
C1 = "#D7301F"
C2 = "#08519C"

with open(SUBTYPE_CLASSIFIER_JSON) as f:
    clf = json.load(f)

ref_genes = clf["genes"]
c1_arr = np.array(clf["centroid_s1"])
c2_arr = np.array(clf["centroid_s2"])


# --- E-MTAB-1980 ---
print("[E-MTAB-1980]")
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
ed1 = cdist(ez[com_em].values, c1_arr[tci_em].reshape(1, -1))[:, 0]
ed2 = cdist(ez[com_em].values, c2_arr[tci_em].reshape(1, -1))[:, 0]
em_labels = np.where(ed1 < ed2, 1, 0)
em_ncc = pd.DataFrame({"sample": list(emat.index), "Subtype": em_labels})

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

em_df = em_ncc.merge(
    clin_em[["sample", "os_m", "event"]], on="sample", how="inner"
).dropna(subset=["os_m", "event"])

em_s1 = em_df[em_df.Subtype == 1]
em_s2 = em_df[em_df.Subtype == 0]
em_lr = logrank_test(em_s1["os_m"], em_s2["os_m"], em_s1["event"], em_s2["event"])
print(
    f"  S1={len(em_s1)} S2={len(em_s2)} events={em_df.event.sum()} "
    f"Log-rank P={em_lr.p_value:.1e}"
)


# --- ICGC RECA-EU ---
print("[ICGC RECA-EU]")
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
id1 = cdist(iez[com_ic].values, c1_arr[tci_ic].reshape(1, -1))[:, 0]
id2 = cdist(iez[com_ic].values, c2_arr[tci_ic].reshape(1, -1))[:, 0]
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

ic_df = pd.DataFrame(ic_rows).dropna()
ic_df["status"] = ic_df["status"].astype(int)

ic_s1 = ic_df[ic_df.Subtype == 1]
ic_s2 = ic_df[ic_df.Subtype == 0]
ic_lr = logrank_test(ic_s1["time"], ic_s2["time"], ic_s1["status"], ic_s2["status"])
print(
    f"  S1={len(ic_s1)} S2={len(ic_s2)} events={ic_df.status.sum()} "
    f"Log-rank P={ic_lr.p_value:.3f}"
)


# --- Figure ---
FS_TITLE = 9
FS_LABEL = 8
FS_TICK = 6.5
FS_LEG = 6.5
FS_PVAL = 7
FS_HR = 6
FS_PANEL = 10

fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.54))
fig.patch.set_facecolor("white")

panels = [
    ("E-MTAB-1980", em_s1, em_s2, em_lr, "os_m", "event"),
    ("ICGC RECA-EU", ic_s1, ic_s2, ic_lr, "time", "status"),
]

for idx, (name, s1, s2, lr, tcol, ecol) in enumerate(panels):
    ax = axes[idx]
    ax.set_facecolor("white")

    ax.text(
        -0.10,
        1.06,
        "AB"[idx],
        transform=ax.transAxes,
        fontsize=FS_PANEL,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax.set_title(name, fontsize=FS_TITLE, fontweight="bold", pad=6)

    kmf = KaplanMeierFitter()
    kmf.fit(s1[tcol], s1[ecol], label=f"S1-like (n={len(s1)})")
    kmf.plot_survival_function(ax=ax, color=C1, lw=1.5)

    kmf.fit(s2[tcol], s2[ecol], label=f"S2-like (n={len(s2)})")
    kmf.plot_survival_function(ax=ax, color=C2, lw=1.5)

    ax.set_xlabel("Time (months)", fontsize=FS_LABEL, fontweight="bold")
    ax.set_ylabel("Overall survival", fontsize=FS_LABEL, fontweight="bold")
    ax.grid(False)
    ax.tick_params(labelsize=FS_TICK)

    if lr.p_value < 0.01:
        p_str = f"P = {lr.p_value:.1e}"
    else:
        p_str = f"P = {lr.p_value:.3f}"

    ax.text(
        0.95,
        0.15,
        p_str,
        transform=ax.transAxes,
        fontsize=FS_PVAL,
        fontweight="bold",
        ha="right",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            alpha=0.9,
            edgecolor="#999",
        ),
    )

    if name == "E-MTAB-1980":
        hr_str = "HR = 4.88 (95% CI: 2.13-11.17)"
    else:
        hr_str = "HR = 1.51 (95% CI: 0.69-3.30)"

    ax.text(
        0.95,
        0.06,
        hr_str,
        transform=ax.transAxes,
        fontsize=FS_HR,
        fontweight="bold",
        ha="right",
        color="#555",
    )

    ax.legend(fontsize=FS_LEG)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout()

out_png = os.path.join(FIG, "External_validation_OS_KM_pubsize.png")
out_svg = os.path.join(FIG, "External_validation_OS_KM_pubsize.svg")
fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
plt.close()
print(f"\nSaved: {os.path.basename(out_png)}")
print(f"       {os.path.basename(out_svg)}")