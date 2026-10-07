# cnv_lipid_remodeling_16genes_check_v2.py
# 16-gene CNA check with BH pooled across the 32 deletion/gain tests,
# matching the 14-gene panel convention. Fisher exact + Spearman CNV-mRNA
# coupling; chromosome arm from the v1 table when available, otherwise
# Ensembl REST (GRCh38).
# Output: results/tables/cnv_lipid_remodeling_16genes_check_v2.csv

import sys, os, io, atexit, json, urllib.request
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests


def _hold():
    if os.environ.get('NOPAUSE'):
        return
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.getch()
    except Exception:
        pass


atexit.register(_hold)

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Expression matrix path: override with KIRC_EXPR_CSV, else look under
# this repo's data/raw/. Adjust to your local layout if needed.
EXPR_CSV = os.environ.get(
    'KIRC_EXPR_CSV',
    os.path.join(PROJ_ROOT, 'data', 'raw', 'KIRC_expr_log2_tpm.csv')
)

OUT = os.path.join(PROJ_ROOT, 'results', 'tables',
                   'cnv_lipid_remodeling_16genes_check_v2.csv')
V1 = os.path.join(PROJ_ROOT, 'results', 'tables',
                  'cnv_lipid_remodeling_16genes_check_v1.csv')
LEGACY = os.path.join(PROJ_ROOT, 'results', 'tables',
                      'cnv_lipid_remodeling_17genes_check_v1.csv')

GENES = {
    "LPCAT1": "tumor-intrinsic", "LPCAT2": "tumor-intrinsic",
    "LPCAT3": "tumor-intrinsic", "LPCAT4": "tumor-intrinsic",
    "MBOAT1": "tumor-intrinsic", "MBOAT2": "tumor-intrinsic",
    "LCLAT1": "tumor-intrinsic",
    "MBOAT7": "tumor-intrinsic", "PLA2G4A": "tumor-intrinsic",
    "PLA2G6": "tumor-intrinsic", "PNPLA8": "tumor-intrinsic",
    "PLA2G4C": "tumor-intrinsic",
    "PLA2G2A": "inflammatory PLA2", "PLA2G4F": "inflammatory PLA2",
    "PLA2G1B": "inflammatory PLA2", "PLA2G2D": "inflammatory PLA2",
}

gistic = pd.read_csv(os.path.join(PROJ_ROOT, "data", "raw",
                                  "TCGA-KIRC.gistic.tsv"),
                     sep="\t", index_col=0)
gistic.columns = [c[:12] for c in gistic.columns]

sub = pd.read_csv(os.path.join(PROJ_ROOT, "data", "processed",
                               "subtype_assignment_balanced.csv"))
sm = dict(zip(sub["Patient"], sub["Subtype"]))

if not os.path.exists(EXPR_CSV):
    raise FileNotFoundError(f"expression matrix not found: {EXPR_CSV}")
expr = pd.read_csv(EXPR_CSV, index_col=0)

common = [c for c in gistic.columns if c in sm and c in expr.columns]
s1 = [c for c in common if sm[c] == "Subtype_1"]
s2 = [c for c in common if sm[c] == "Subtype_2"]
print(f"matched patients: {len(common)}  (S1={len(s1)}, S2={len(s2)})")

CEN = {
    "1": 123.4, "2": 93.3, "3": 90.9, "4": 50.3, "5": 48.6, "6": 61.0,
    "7": 60.0, "8": 45.6, "9": 43.0, "10": 40.2, "11": 53.7, "12": 36.6,
    "13": 18.0, "14": 17.6, "15": 19.0, "16": 37.2, "17": 25.1, "18": 18.5,
    "19": 26.2, "20": 28.0, "21": 11.4, "22": 14.7, "X": 61.0, "Y": 11.0,
}


def arm_of(gene):
    try:
        url = (f"https://rest.ensembl.org/lookup/human/{gene}"
               f"?content-type=application/json")
        req = urllib.request.Request(url, headers={"User-Agent": "check"})
        with urllib.request.urlopen(req, timeout=10) as r:
            d = json.loads(r.read().decode())
        chrom = str(d["seq_region_name"])
        mid = (d["start"] + d["end"]) / 2 / 1e6
        if chrom not in CEN:
            return "NA"
        return f"{chrom}{'p' if mid < CEN[chrom] else 'q'}"
    except Exception:
        return "NA"


arms = {g: "NA" for g in GENES}

legacy_arms = {}
if os.path.exists(LEGACY):
    leg = pd.read_csv(LEGACY)
    legacy_arms = dict(zip(leg["gene"], leg["arm"]))

if os.path.exists(V1):
    v1 = pd.read_csv(V1)
    if "arm" in v1.columns:
        for g, a in zip(v1["gene"], v1["arm"]):
            if (g in arms and isinstance(a, str) and a
                    and not a.startswith("NA")):
                arms[g] = a

for g, a in arms.items():
    if a == "NA" and g in legacy_arms:
        arms[g] = legacy_arms[g]
for g, a in arms.items():
    if a == "NA":
        arms[g] = arm_of(g)

rows_del, rows_amp, rows_core = [], [], []
for g, group in GENES.items():
    if g not in gistic.index:
        rows_del.append({"gene": g})
        rows_amp.append({"gene": g})
        continue

    v = gistic.loc[g, common].astype(float)
    e = expr.loc[g, common].astype(float) if g in expr.index else None

    for kind, thr, dst in (("del", -1, rows_del), ("amp", 1, rows_amp)):
        if kind == "del":
            alt1 = int((v[s1] <= thr).sum())
            alt2 = int((v[s2] <= thr).sum())
        else:
            alt1 = int((v[s1] >= thr).sum())
            alt2 = int((v[s2] >= thr).sum())

        t = [[alt1, len(s1) - alt1], [alt2, len(s2) - alt2]]
        p = fisher_exact(t, alternative="two-sided")[1]
        dst.append({
            "gene": g,
            f"S1_{kind}_pct": round(100 * alt1 / len(s1), 1),
            f"S2_{kind}_pct": round(100 * alt2 / len(s2), 1),
            f"p_{kind}": p,
        })

    if e is not None:
        rho, pr = spearmanr(v, e, nan_policy="omit")
    else:
        rho, pr = np.nan, np.nan

    rows_core.append({
        "gene": g, "group": group, "arm": arms[g],
        "mRNA_S1_mean": round(e[s1].mean(), 3) if e is not None else np.nan,
        "mRNA_S2_mean": round(e[s2].mean(), 3) if e is not None else np.nan,
        "spearman_rho_CNV_mRNA":
            round(float(rho), 3) if np.isfinite(rho) else np.nan,
        "rho_p": float(pr) if np.isfinite(pr) else np.nan,
    })

res = pd.DataFrame(rows_core)
if res.empty:
    raise RuntimeError("no genes processed; check GISTIC index")

d = pd.DataFrame(rows_del)
a = pd.DataFrame(rows_amp)

res = res.merge(
    d.reindex(columns=["gene", "S1_del_pct", "S2_del_pct", "p_del"]),
    on="gene", how="left",
).merge(
    a.reindex(columns=["gene", "S1_amp_pct", "S2_amp_pct", "p_amp"]),
    on="gene", how="left",
)

# BH pooled over all available del/amp tests (single family, matches the
# 14-gene panel convention: 28 tests for 14 genes x 2 directions).
mask = res["p_del"].notna() & res["p_amp"].notna()
if mask.any():
    pool = multipletests(
        np.concatenate([res.loc[mask, "p_del"].values,
                        res.loc[mask, "p_amp"].values]),
        method="fdr_bh")[1]
    k = int(mask.sum())
    res.loc[mask, "p_del_adj"] = pool[:k]
    res.loc[mask, "p_amp_adj"] = pool[k:]
else:
    res["p_del_adj"] = np.nan
    res["p_amp_adj"] = np.nan

res["rho_p_adj"] = np.nan
mask = res["rho_p"].notna()
if mask.any():
    res.loc[mask, "rho_p_adj"] = multipletests(
        res.loc[mask, "rho_p"].astype(float), method="fdr_bh")[1]

res["del_sig"] = np.where(res["p_del_adj"] < 0.05, "*", "")
res["amp_sig"] = np.where(res["p_amp_adj"] < 0.05, "*", "")
res["cnv_mRNA_coupled"] = np.where(res["rho_p_adj"] < 0.05, "yes", "no")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
res.to_csv(OUT, index=False)

print(res.round(4).to_string(index=False))
print(f"\nSaved: {OUT}")