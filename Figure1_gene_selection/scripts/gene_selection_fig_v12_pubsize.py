#!/usr/bin/env python3
"""Generate Figure 1: gene selection overview."""

import io
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import norm

warnings.filterwarnings("ignore")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.linewidth"] = 0.8

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
WORK = PROJECT / "data"
FIG = PROJECT / "results"
FIG.mkdir(parents=True, exist_ok=True)

with open(WORK / "msigdb_lipid_metadata.json") as f:
    meta = json.load(f)
with open(WORK / "go_dedup_stats.json") as f:
    go_stats = json.load(f)
with open(WORK / "cox_screening" / "cox_summary.json") as f:
    cox_stats = json.load(f)

strict = pd.read_csv(WORK / "gene_selection" / "selected_genes_strict.csv")
cox_all = pd.read_csv(WORK / "cox_screening" / "cox_results_all.csv")

cox_all["p_value"] = pd.to_numeric(
    cox_all["p_value"].astype(str).str.replace(r"[\[\]]", "", regex=True),
    errors="coerce",
)
cox_all["coef"] = pd.to_numeric(cox_all["coef"], errors="coerce")
cox_all["FDR"] = pd.to_numeric(
    cox_all["FDR"].astype(str).str.replace(r"[\[\]]", "", regex=True),
    errors="coerce",
)
if "se" in cox_all.columns:
    cox_all["se"] = pd.to_numeric(
        cox_all["se"].astype(str).str.replace(r"[\[\]]", "", regex=True),
        errors="coerce",
    )

with open(WORK / "consensus_cluster" / "balanced_genes.txt") as f:
    balanced = [line.strip() for line in f if line.strip()]

bal_df = strict[strict["gene"].isin(balanced)]
layer_counts = bal_df["functional_layer"].value_counts()

top_cox = cox_all.nsmallest(15, "p_value").copy()
if "se" in top_cox.columns and top_cox["se"].notna().any():
    top_cox["SE"] = top_cox["se"]
else:
    top_cox["SE"] = abs(top_cox["coef"]) / norm.ppf(
        1 - top_cox["p_value"] / 2
    ).clip(1e-10)

top_cox["HR"] = np.exp(top_cox["coef"])
top_cox["CI_low"] = np.exp(top_cox["coef"] - 1.96 * top_cox["SE"])
top_cox["CI_high"] = np.exp(top_cox["coef"] + 1.96 * top_cox["SE"])

FS_PANEL = 10
FS_LABEL = 8
FS_TICK = 7
FS_GENE = 7
FS_LEG = 7
FS_CNT = 7
FS_VOLC = 7

fig = plt.figure(figsize=(7.09, 4.4))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.55, 1.0],
    hspace=0.42,
    wspace=0.22,
    left=0.09,
    right=0.975,
    top=0.95,
    bottom=0.09,
)

ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_facecolor("white")
log_p = -np.log10(cox_all["FDR"].clip(1e-300))
sig_m = cox_all["FDR"] < 0.05
n_red = int(sig_m.sum())
n_gray = int((~sig_m).sum())

ax_a.scatter(cox_all["coef"][~sig_m], log_p[~sig_m], c="#BDC3C7", s=2, alpha=0.25)
ax_a.scatter(
    cox_all["coef"][sig_m],
    log_p[sig_m],
    c="#E74C3C",
    s=6,
    alpha=0.6,
    edgecolors="none",
)
ax_a.text(
    0.95, 0.95,
    f"Red: significant genes (n={n_red})\nGrey: non-significant genes (n={n_gray})",
    transform=ax_a.transAxes,
    fontsize=FS_VOLC,
    va="top",
    ha="right",
    color="#333",
)
ax_a.axhline(y=-np.log10(0.05), color="#333", ls="--", lw=0.5, alpha=0.3)
ax_a.set_xlabel("log(HR)", fontsize=FS_LABEL, fontweight="bold")
ax_a.set_ylabel(r"$-\mathregular{log_{10}}$ (FDR)", fontsize=FS_LABEL, fontweight="bold")
ax_a.tick_params(labelsize=FS_TICK)
ax_a.text(-0.12, 1.04, "A", transform=ax_a.transAxes, fontsize=FS_PANEL, fontweight="bold", va="bottom")

ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_facecolor("white")
genes_d = list(top_cox["gene"])[::-1]
hrs = list(top_cox["HR"])[::-1]
lows = list(top_cox["CI_low"])[::-1]
highs = list(top_cox["CI_high"])[::-1]

for i, (hr, lo, hi) in enumerate(zip(hrs, lows, highs)):
    color = "#B71C1C" if hr > 1 else "#0D47A1"
    ax_b.errorbar(hr, i, xerr=[[hr - lo], [hi - hr]], fmt="o", color=color, capsize=2, ms=5, lw=0.8)

ax_b.axvline(x=1, color="#333", ls="--", lw=0.8)
legend_b = [
    Patch(facecolor="#B71C1C", edgecolor="none", label="Risk (HR > 1)"),
    Patch(facecolor="#0D47A1", edgecolor="none", label="Protective (HR < 1)"),
]
ax_b.legend(
    handles=legend_b,
    loc="lower right",
    bbox_to_anchor=(0.98, 0.02),
    frameon=True,
    framealpha=0.95,
    fontsize=FS_LEG,
)
ax_b.set_xlim(0.4, 1.8)
ax_b.set_yticks(range(len(genes_d)))
ax_b.set_yticklabels(genes_d, fontsize=FS_GENE, fontweight="bold", fontstyle="italic")
ax_b.set_xlabel("HR (95% CI)", fontsize=FS_LABEL, fontweight="bold")
ax_b.tick_params(labelsize=FS_TICK)
ax_b.text(-0.12, 1.04, "B", transform=ax_b.transAxes, fontsize=FS_PANEL, fontweight="bold", va="bottom")

ax_c = fig.add_subplot(gs[1, :])
ax_c.set_facecolor("white")

fig.canvas.draw()
x_left_c = ax_a.transData.transform((0.1, 0))[0] / (fig.dpi * fig.get_figwidth())
x_right_c = ax_b.get_position().x1
pos_c = ax_c.get_position()
ax_c.set_position([x_left_c, pos_c.y0, x_right_c - x_left_c, pos_c.height])

short_names = {
    "Core_Lipid_Enzyme": "Core lipid enzymes",
    "Signaling_Lipoprotein": "Signaling / Lipoprotein-associated genes",
    "Transcriptional_Regulator": "Transcriptional regulators",
    "Other": "Other lipid metabolism genes",
    "Lipid_Droplet_Storage": "Lipid droplet storage factors",
    "Ferroptosis_ROS": "Ferroptosis / ROS-related genes",
}

layers = layer_counts.index.tolist()[::-1]
vals_e = layer_counts.values.tolist()[::-1]
labels_e = [short_names.get(layer, layer) for layer in layers]
bar_colors = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]
bar_colors_ordered = bar_colors[:len(layers)][::-1]

ax_c.barh(
    range(len(layers)),
    vals_e,
    color=bar_colors_ordered,
    alpha=0.85,
    edgecolor="white",
)
ax_c.set_yticks(range(len(layers)))
ax_c.set_yticklabels(labels_e, fontsize=FS_GENE, fontweight="bold")
ax_c.set_xlabel("Number of genes", fontsize=FS_LABEL, fontweight="bold")
ax_c.tick_params(labelsize=FS_TICK)
ax_c.set_xlim(0, max(vals_e) * 1.15)
ax_c.text(-0.12, 1.04, "C", transform=ax_c.transAxes, fontsize=FS_PANEL, fontweight="bold", va="bottom")

for i, v in enumerate(vals_e):
    ax_c.text(v + 0.3, i, str(v), va="center", fontsize=FS_CNT, fontweight="bold", color=bar_colors_ordered[i])

fig.savefig(
    FIG / "Figure_1_gene_selection_v12_pubsize.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    pad_inches=0.05,
)
fig.savefig(
    FIG / "Figure_1_gene_selection_v12_pubsize.svg",
    format="svg",
    bbox_inches="tight",
    facecolor="white",
    pad_inches=0.05,
)
plt.close()
print("saved Figure_1_gene_selection_v12_pubsize.png + SVG")