# fig_driver_3p_del_FAO9_v1.py
# Two-panel supplementary figure:
#   A: 3p driver deletion rates, S1 vs S2 (GISTIC gene-level)
#   B: FAO9 score by driver mutation status
# Output: results/figures/Supplementary_Fig_driver_3p_del_FAO9_v1.{png,svg}

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig4_params as P

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']

DEL_CLR = '#B71C1C'
WT_CLR = '#D9D9D9'
MUT_CLR = '#B71C1C'
SPINE_LW = 0.8
FS_TICK, FS_LABEL, FS_STAR, FS_P = 10, 11, 10, 8
BAR_H = 0.22
GENE_SPAN = 1.30
BOX_W = 0.30
X_OFF = 0.19

GENES = P.MUT_GENES_ORDER


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    return np.minimum(adj, 1.0)[np.argsort(order)]


def stars(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


# ---- Panel A ----
del_df = pd.read_csv(os.path.join(
    P.PROJ_ROOT, 'results', 'tables', 'driver_gistic_del_v1.csv'))
del_genes = del_df['gene'].tolist()
s1_del = del_df['S1_del_pct'].to_numpy()
s2_del = del_df['S2_del_pct'].to_numpy()
del_fdr = bh_adjust(del_df['p_fisher_del_S1vsS2'].to_numpy())
n1 = int(del_df['S1_n'].iloc[0])
n2 = int(del_df['S2_n'].iloc[0])

# ---- Panel B ----
fao9 = pd.read_csv(os.path.join(P.PROC_DIR, 'fao9_immune_input_tcga_v1.csv'))
fao9['Patient'] = fao9['sample'].str[:12]
fao9_map = dict(zip(fao9['Patient'], fao9['FAO9']))

maf = pd.read_csv(os.path.join(P.RAW_DIR, 'TCGA-KIRC.merged.maf'),
                  sep='\t', comment='#', low_memory=False)
maf = maf[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification']].copy()
maf = maf[maf['Variant_Classification'].isin(P.NONSYN)].copy()
maf['Patient'] = maf['Tumor_Sample_Barcode'].str[:12]
cohort = sorted(set(maf['Patient']) & set(fao9_map))

mut_tab = pd.read_csv(os.path.join(
    P.PROJ_ROOT, 'results', 'tables', 'driver_mut_fao_v1.csv')).set_index('gene')

values = {}
for gene in GENES:
    mut = set(maf.loc[maf['Hugo_Symbol'] == gene, 'Patient'])
    v_mut = [fao9_map[p] for p in cohort if p in mut]
    v_wt = [fao9_map[p] for p in cohort if p not in mut]
    _, p_mwu = mannwhitneyu(v_mut, v_wt, alternative='two-sided')
    assert abs(p_mwu - float(mut_tab.loc[gene, 'p_mwu_FAO9'])) < 1e-9, \
        f'{gene}: MWU P mismatch vs driver_mut_fao_v1.csv'
    values[gene] = (np.asarray(v_wt), np.asarray(v_mut))

print('MWU P values cross-checked against driver_mut_fao_v1.csv: all match')

# ---- figure ----
fig = plt.figure(figsize=(7.087, 5.6), dpi=600)
fig.patch.set_facecolor('white')
gs = fig.add_gridspec(2, 1, height_ratios=[2.1, 2.9])
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])
for ax in (axA, axB):
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)

# Panel A
n_genes = len(del_genes)
y_gene = (n_genes - 1 - np.arange(n_genes)) * GENE_SPAN
y_off = np.array([0.15, -0.15])

axA.barh(y_gene + y_off[0], s1_del, height=BAR_H,
         color=DEL_CLR, alpha=0.45, zorder=3)
axA.barh(y_gene + y_off[1], s2_del, height=BAR_H,
         color=DEL_CLR, alpha=0.85, zorder=3)

xmax = max(s1_del.max(), s2_del.max()) * 1.25
for i, f in enumerate(del_fdr):
    s = stars(f)
    if s:
        axA.text(s1_del[i] + xmax * 0.015, y_gene[i] + y_off[0], s,
                 va='center', ha='left', fontsize=FS_STAR,
                 fontweight='bold', color=DEL_CLR)
        axA.text(s2_del[i] + xmax * 0.015, y_gene[i] + y_off[1], s,
                 va='center', ha='left', fontsize=FS_STAR,
                 fontweight='bold', color=DEL_CLR)

axA.set_yticks(y_gene)
axA.set_yticklabels(del_genes, fontsize=FS_TICK,
                    fontstyle='italic', fontweight='bold')
for tl in axA.get_yticklabels():
    tl.set_color('black')
axA.set_xlabel('Patients with CNA (%)', fontsize=FS_LABEL, fontweight='bold')
axA.set_xlim(0, xmax)
axA.set_ylim(y_gene[-1] - 0.70, y_gene[0] + 0.70)
axA.tick_params(axis='both', labelsize=FS_TICK)
axA.grid(True, alpha=0.10, axis='x')

axA.legend(handles=[
    Patch(facecolor=DEL_CLR, alpha=0.45, label=f'S1 deletion (n = {n1})'),
    Patch(facecolor=DEL_CLR, alpha=0.85, label=f'S2 deletion (n = {n2})'),
], fontsize=7.5, loc='lower left', frameon=False,
   bbox_to_anchor=(0.0, 1.02), ncol=2,
   handlelength=1.0, handleheight=0.8, columnspacing=0.8, handletextpad=0.5)

# Panel B
rng = np.random.default_rng(0)
all_vals = np.concatenate([np.concatenate(v) for v in values.values()])
vmin, vmax = all_vals.min(), all_vals.max()
vspan = vmax - vmin
axB.set_ylim(vmin - 0.10 * vspan, vmax + 0.22 * vspan)

for i, gene in enumerate(GENES):
    v_wt, v_mut = values[gene]

    for x, v, face, edge in [(i - X_OFF, v_wt, WT_CLR, '0.45'),
                             (i + X_OFF, v_mut, MUT_CLR, MUT_CLR)]:
        bp = axB.boxplot(v, positions=[x], widths=BOX_W, patch_artist=True,
                         showfliers=False, zorder=3,
                         medianprops=dict(color='black', linewidth=1.2),
                         whiskerprops=dict(color='0.30', linewidth=0.8),
                         capprops=dict(color='0.30', linewidth=0.8))
        bp['boxes'][0].set(facecolor=face, edgecolor=edge, linewidth=0.8)

    jw = rng.uniform(-0.085, 0.085, len(v_wt))
    axB.scatter(i - X_OFF + jw, v_wt, s=3.5, color='0.35',
                alpha=0.22, linewidths=0, zorder=2)
    jm = rng.uniform(-0.085, 0.085, len(v_mut))
    axB.scatter(i + X_OFF + jm, v_mut, s=8, color=MUT_CLR,
                alpha=0.55, linewidths=0, zorder=4)

    p_bh = float(mut_tab.loc[gene, 'p_bh_mwu'])
    axB.text(i, vmax + 0.13 * vspan, f'P = {p_bh:.2f} (ns)',
             ha='center', va='center', fontsize=FS_P, color='black')
    axB.text(i, vmin - 0.045 * vspan, f'{len(v_wt)} / {len(v_mut)}',
             ha='center', va='center', fontsize=6.5, color='0.35')

axB.set_xticks(range(len(GENES)))
axB.set_xticklabels(GENES, fontsize=FS_TICK,
                    fontstyle='italic', fontweight='bold')
axB.tick_params(axis='both', labelsize=FS_TICK, length=3)
axB.set_ylabel('FAO9 score', fontsize=FS_LABEL, fontweight='bold')
axB.grid(True, axis='y', alpha=0.10)

axB.legend(handles=[
    Patch(facecolor=WT_CLR, edgecolor='0.45', label='Wild-type'),
    Patch(facecolor=MUT_CLR, edgecolor=MUT_CLR, label='Mutant'),
], fontsize=7.5, loc='lower left', frameon=False,
   bbox_to_anchor=(0.0, 1.02), ncol=2,
   handlelength=1.0, handleheight=0.8, columnspacing=0.8, handletextpad=0.5)

fig.subplots_adjust(left=0.115, right=0.965, top=0.905,
                    bottom=0.085, hspace=0.62)

for ax, letter in ((axA, 'A'), (axB, 'B')):
    bb = ax.get_position()
    fig.text(bb.x0 - 0.045, bb.y1 + 0.035, letter, fontsize=13,
             fontweight='bold', va='top', ha='left')

out_png = os.path.join(P.OUT_DIR, 'Supplementary_Fig_driver_3p_del_FAO9_v1.png')
out_svg = os.path.join(P.OUT_DIR, 'Supplementary_Fig_driver_3p_del_FAO9_v1.svg')
fig.savefig(out_png, dpi=600)
fig.savefig(out_svg)

print('saved:', out_png)
print('saved:', out_svg)
print('Panel A stars:', {g: stars(f) for g, f in zip(del_genes, del_fdr)})
import os
print("\n" + "=" * 60)
print()
print("=" * 60)
if os.name == "nt":
    import msvcrt
    msvcrt.getch()
else:
    input()