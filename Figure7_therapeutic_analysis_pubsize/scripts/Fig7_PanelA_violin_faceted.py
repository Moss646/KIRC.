#!/usr/bin/env python3
# Fig7_PanelA_violin_faceted.py
# 3x3 faceted violin/box plot of predicted ln(IC50) per drug (S1 vs S2),
# with Mann-Whitney U + BH-FDR significance stars.
# Output: results/figures/Fig7_PanelA_violin_faceted.{png,svg}

import sys, os, io, atexit
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu


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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(REPO, 'data', 'processed')
OUTDIR = os.path.join(REPO, 'results', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.unicode_minus": False,
})

C1, C2 = '#C65A5A', '#4C78A8'
np.random.seed(42)

ic50_raw = pd.read_csv(os.path.join(PROC, 'oncoPredict_calcPheno_ic50.csv'),
                       index_col=0)
sub = pd.read_csv(os.path.join(PROC, 'subtype_assignment_balanced.csv'))

sub_map = dict(zip(sub['Patient'].str.replace('-', '.').str[:15],
                   sub['Subtype']))


def norm_id(x):
    x = str(x).replace('-', '.').replace('_', '.')
    if not x.startswith('TCGA'):
        x = 'TCGA.' + x
    return x[:15]


ic50_raw['Group'] = ic50_raw.index.map(lambda x: sub_map.get(norm_id(x)))
ic50_raw = ic50_raw.dropna(subset=['Group']).copy()
ic50_raw['Group'] = ic50_raw['Group'].replace({'Subtype_1': 'S1',
                                               'Subtype_2': 'S2'})

drugs9 = ['Axitinib', 'Sorafenib', 'Rapamycin', 'Cediranib', 'Crizotinib',
          'Dasatinib', 'Trametinib', 'AZD2014', 'AZD8055']

S1 = ic50_raw[ic50_raw.Group == 'S1']
S2 = ic50_raw[ic50_raw.Group == 'S2']
print(f'Samples: S1={len(S1)}, S2={len(S2)}, total={len(ic50_raw)}')

pvals = []
for d in drugs9:
    _, p = mannwhitneyu(S1[d].values, S2[d].values, alternative='two-sided')
    pvals.append(p)

pvals_arr = np.array(pvals)
order = np.argsort(pvals_arr)
ranked = pvals_arr[order]
m = len(ranked)
fdr_ranked = ranked * m / (np.arange(m) + 1)
fdr_ranked = np.minimum.accumulate(fdr_ranked[::-1])[::-1]
fdr_ranked = np.clip(fdr_ranked, 0, 1)
fdr_full = np.empty_like(fdr_ranked)
fdr_full[order] = fdr_ranked
fdr_map = dict(zip(drugs9, fdr_full))


def star(fdr):
    if fdr < 0.001:
        return '***'
    if fdr < 0.01:
        return '**'
    if fdr < 0.05:
        return '*'
    return 'ns'


print('\nDrug        MWU p        BH-FDR      star')
for d, p, f in zip(drugs9, pvals, fdr_full):
    print(f'{d:<11} {p:.2e}   {f:.2e}   {star(f)}')

fig, axes = plt.subplots(3, 3, figsize=(14, 14))
fig.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.06,
                    wspace=0.18, hspace=0.45)

for idx, drug in enumerate(drugs9):
    ax = axes[idx // 3, idx % 3]
    data_s1 = S1[drug].values
    data_s2 = S2[drug].values

    positions = [1, 2]
    data_list = [data_s1, data_s2]
    colors = [C1, C2]

    parts = ax.violinplot(data_list, positions=positions, showmeans=False,
                          showmedians=False, showextrema=False, widths=0.7)
    for pc, c in zip(parts['bodies'], colors):
        pc.set_facecolor(c)
        pc.set_alpha(0.25)
        pc.set_edgecolor(c)
        pc.set_linewidth(1.5)

    bp = ax.boxplot(data_list, positions=positions, widths=0.25,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='#333333', lw=2.0),
                    whiskerprops=dict(lw=1.5, color='#333333'),
                    capprops=dict(lw=1.5, color='#333333'),
                    boxprops=dict(lw=1.5, color='#333333'))
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)

    for pos, vals, c in zip(positions, data_list, colors):
        jitter = np.random.normal(0, 0.06, len(vals))
        ax.scatter(np.ones(len(vals)) * pos + jitter, vals, alpha=0.15, s=8,
                   color=c, edgecolors='none', zorder=3)

    ymax = max(data_s1.max(), data_s2.max())
    ymin = min(data_s1.min(), data_s2.min())
    bracket_y = ymax + (ymax - ymin) * 0.12
    ax.plot([1, 1, 2, 2],
            [bracket_y - 0.01, bracket_y, bracket_y, bracket_y - 0.01],
            color='#333333', lw=1.2)
    s = star(fdr_map[drug])
    fs = 14 if s != 'ns' else 11
    ax.text(1.5, bracket_y + (ymax - ymin) * 0.02, s, ha='center',
            va='bottom', fontsize=fs, fontweight='bold', color='#333333')

    ax.set_xticks(positions)
    ax.set_xticklabels(['S1', 'S2'], fontsize=12, fontweight='bold')
    ax.set_title(drug, fontsize=14, fontweight='bold', pad=8)
    ax.tick_params(axis='y', labelsize=11)

    y_pad = (ymax - ymin) * 0.05
    ax.set_ylim(ymin - y_pad, bracket_y + (ymax - ymin) * 0.10)

    for sp in ['top', 'bottom', 'left', 'right']:
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_color('#444444')
        ax.spines[sp].set_linewidth(1.5)

fig.text(0.025, 0.5, 'Predicted ln(IC$_{50}$)', va='center',
         rotation='vertical', fontsize=14, fontweight='bold')

legend_handles = [Patch(facecolor=C1, alpha=0.6,
                        label='S1 (n={})'.format(len(S1))),
                  Patch(facecolor=C2, alpha=0.6,
                        label='S2 (n={})'.format(len(S2)))]
fig.legend(handles=legend_handles, loc='upper right', fontsize=11,
           framealpha=0.9, bbox_to_anchor=(0.97, 0.99))

fig.text(0.5, 0.02,
         'Mann\u2013Whitney U test with BH-FDR correction:  '
         '* FDR < 0.05,  ** FDR < 0.01,  *** FDR < 0.001',
         ha='center', fontsize=10, color='#555555')

out_png = os.path.join(OUTDIR, 'Fig7_PanelA_violin_faceted.png')
out_svg = os.path.join(OUTDIR, 'Fig7_PanelA_violin_faceted.svg')
fig.savefig(out_png, dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
fig.savefig(out_svg, dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none',
            metadata={'Date': '2026-08-12T00:00:00'})
plt.close(fig)

print(f'\nSaved: {out_png}')
print(f'Saved: {out_svg}')