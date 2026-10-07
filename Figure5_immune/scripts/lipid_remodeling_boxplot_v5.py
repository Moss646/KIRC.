# lipid_remodeling_boxplot_v5.py
# S1 vs S2 expression of 16 lipid-remodeling genes in two functional groups.
# Boxes show log2(TPM+1); stars show limma-voom FDR from the Figure3 table.
# Output: results/tables/lipid_remodeling_16genes_stats_v5.csv
#         results/figures/lipid_remodeling_16genes_S1S2_v5.{png,svg}

import sys, os, io, atexit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _hold():
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_ROOT = os.environ.get('MEL_DATA_ROOT', PROJ_ROOT)

EXPR_CSV = os.path.join(DATA_ROOT, 'data', 'tcga_kirc',
                        'KIRC_expr_log2_tpm.csv')
SUB_CSV = os.path.join(DATA_ROOT, 'data', 'processed',
                       'subtype_assignment_balanced.csv')
LIMMA_CSV = os.environ.get(
    'LIMMA_CSV',
    os.path.join(DATA_ROOT, 'data', 'processed',
                 'deg_limma_voom_with_symbols.csv')
)
TAB_OUT = os.path.join(PROJ_ROOT, 'results', 'tables',
                       'lipid_remodeling_16genes_stats_v5.csv')
FIG_OUT = os.path.join(PROJ_ROOT, 'results', 'figures',
                       'lipid_remodeling_16genes_S1S2_v5')

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

GROUPS = [
    ("Tumor-intrinsic",
     ["LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4",
      "MBOAT1", "MBOAT2", "LCLAT1", "MBOAT7",
      "PLA2G4A", "PLA2G6", "PNPLA8", "PLA2G4C"]),
    ("Inflammatory PLA2",
     ["PLA2G2A", "PLA2G4F", "PLA2G1B", "PLA2G2D"]),
]
GENES = [g for _, gs in GROUPS for g in gs]
DISPLAY = {"LCLAT1": "MBOAT5\n(LCLAT1)"}

expr = pd.read_csv(EXPR_CSV, index_col=0)
sub = pd.read_csv(SUB_CSV)
sm = sub.set_index("Patient")["Subtype"]
common = [c for c in expr.columns if c in sm.index]
s1_ids = [c for c in common if sm[c] == "Subtype_1"]
s2_ids = [c for c in common if sm[c] == "Subtype_2"]

missing = [g for g in GENES if g not in expr.index]
if missing:
    raise SystemExit(f"genes missing from matrix: {missing}")

limma = pd.read_csv(LIMMA_CSV)
limma_idx = limma.set_index("gene_symbol")
missing_limma = [g for g in GENES if g not in limma_idx.index]
if missing_limma:
    raise SystemExit(f"genes missing from limma table: {missing_limma}")


def fdr_to_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


rows = []
for gname, genes in GROUPS:
    for g in genes:
        v1 = expr.loc[g, s1_ids].astype(float)
        v2 = expr.loc[g, s2_ids].astype(float)
        lr = limma_idx.loc[g]
        fdr = float(lr["FDR"])
        rows.append({
            "gene": g,
            "group": gname,
            "S1_mean": v1.mean(),
            "S2_mean": v2.mean(),
            "diff_S2_minus_S1": v2.mean() - v1.mean(),
            "limma_logFC": float(lr["logFC"]),
            "limma_FDR": fdr,
            "limma_category": lr["category"],
            "star": fdr_to_stars(fdr),
        })

res = pd.DataFrame(rows)
os.makedirs(os.path.dirname(TAB_OUT), exist_ok=True)
res.to_csv(TAB_OUT, index=False)

C1, C2 = "#C65A5A", "#4C78A8"
fig, ax = plt.subplots(figsize=(7.09, 3.4))
rng = np.random.default_rng(1)

n = len(GENES)
sep_x = None
for k, g in enumerate(GENES):
    if k == len(GROUPS[0][1]):
        sep_x = k - 0.5
    v1 = expr.loc[g, s1_ids].astype(float).values
    v2 = expr.loc[g, s2_ids].astype(float).values
    for off, vals, col in ((-0.18, v1, C1), (0.18, v2, C2)):
        ax.boxplot([vals], positions=[k + off], widths=0.30,
                   patch_artist=True, showfliers=False,
                   medianprops=dict(color="black", lw=1.0),
                   boxprops=dict(facecolor=col, alpha=0.55, lw=0.8),
                   whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8))
        jitter = rng.normal(k + off, 0.05, size=len(vals))
        ax.scatter(jitter, vals, s=2.5, color=col, alpha=0.45, lw=0)
    st = res.loc[res["gene"] == g, "star"].iloc[0]
    if st:
        ymax = max(v1.max(), v2.max())
        ymin = min(v1.min(), v2.min())
        ax.text(k, ymax + 0.12 * (ymax - ymin), st, ha="center", va="bottom",
                fontsize=8, fontweight="bold", color="red")

ax.axvline(x=sep_x, color="#999999", lw=1.2, ls="--", alpha=0.6, zorder=0)

ax.text((n - 4 - 1) / 2, 1.045, GROUPS[0][0],
        transform=ax.get_xaxis_transform(),
        ha="center", va="bottom", fontsize=8.5, color="#555555")
ax.text(n - 2, 1.045, GROUPS[1][0], transform=ax.get_xaxis_transform(),
        ha="center", va="bottom", fontsize=8.5, color="#555555")

ax.set_xticks(range(n))
ax.set_xticklabels([DISPLAY.get(g, g) for g in GENES], rotation=40,
                   ha="right", fontsize=8)
for t in ax.get_xticklabels():
    t.set_style("italic")
ax.set_ylabel("Expression, log2(TPM+1)", fontsize=9)
ax.set_xlim(-0.6, n - 0.4)
for sp in ax.spines.values():
    sp.set_linewidth(0.8)

h = [plt.Rectangle((0, 0), 1, 1, fc=C1, alpha=0.55),
     plt.Rectangle((0, 0), 1, 1, fc=C2, alpha=0.55)]
ax.legend(h, [f"S1 (n={len(s1_ids)})", f"S2 (n={len(s2_ids)})"],
          loc="upper left", frameon=False, fontsize=8.5)

os.makedirs(os.path.dirname(FIG_OUT), exist_ok=True)
fig.savefig(FIG_OUT + ".png", dpi=600, bbox_inches="tight")
fig.savefig(FIG_OUT + ".svg", bbox_inches="tight")
print("saved:", TAB_OUT)
print("saved:", FIG_OUT + ".png / .svg")
print(res.round(4).to_string(index=False))