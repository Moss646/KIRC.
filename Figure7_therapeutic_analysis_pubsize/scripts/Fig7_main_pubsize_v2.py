# Fig7_main_pubsize_v2.py
# Figure 7 publication-size version v2 (180 x 210 mm @ 600 dpi).
# Panels: A drug sensitivity | B Cohen's d | C residual Cohen's d
#         D responder rate  | E TIDE scores | F immune markers
# Output: results/figures/Figure_7_therapeutic_analysis_pubsize_v2.{png,svg}

import sys, os, io, atexit
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')


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

FS_TITLE = 9
FS_LABEL = 8
FS_TICK = 7
FS_TICK_SM = 6.5
FS_TICK_BC = 6.0
FS_LEG = 6.5
FS_VAL = 6.5
FS_STAR = 8
FS_STAR_A = 6.5
FS_STAR_D = 9
LW = 0.8

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "svg.fonttype": "none",
    "axes.linewidth": LW,
    "svg.hashsalt": "fig7pubsize_v2",
})

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_ROOT = os.path.join(PROJ_ROOT, 'data')

PROC = os.path.join(DATA_ROOT, 'processed')
FIG_OUT = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(FIG_OUT, exist_ok=True)

C1, C2 = '#C65A5A', '#4C78A8'
np.random.seed(42)


def star(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


tide = pd.read_csv(os.path.join(PROC, 'TIDE_official_results.csv'), index_col=0)
tide.index = tide.index.str.strip()
for c in tide.select_dtypes(include='bool').columns:
    tide[c] = tide[c].astype(float)

sub = pd.read_csv(os.path.join(PROC, 'subtype_assignment_balanced.csv'))
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
tide['Group'] = tide.index.map(
    lambda x: sub_map.get(x, sub_map.get(x[:12], None)))
tide = tide.dropna(subset=['Group'])
tide['Group'] = tide['Group'].replace({'Subtype_1': 'S1', 'Subtype_2': 'S2'})

s1_n = (tide.Group == 'S1').sum()
s2_n = (tide.Group == 'S2').sum()
s1_resp = tide[tide.Group == 'S1']['Responder'].astype(bool).sum()
s2_resp = tide[tide.Group == 'S2']['Responder'].astype(bool).sum()

# Panel A: official oncoPredict::calcPhenotype values. Drug order is
# determined downstream by the residual Cohen's d (ABC_ORDER).
onco_off = pd.read_csv(os.path.join(PROC, 'oncoPredict_calcPheno_results.csv'))
onco_off['delta_lnIC50'] = onco_off['S1_mean'] - onco_off['S2_mean']

drug_names = onco_off['Drug'].tolist()
delta_display = onco_off['delta_lnIC50'].tolist()
p_vals = onco_off['p'].tolist()
fdr_vals = onco_off['FDR'].tolist()

# Per-sample IC50 for panels B and C
ic50_raw = pd.read_csv(os.path.join(PROC, 'oncoPredict_calcPheno_ic50.csv'),
                       index_col=0)


def norm_id(s):
    return str(s).replace('.', '-')


sub_map2 = dict(zip(sub['Patient'].map(norm_id), sub['Subtype']))
ic50_raw['G'] = ic50_raw.index.map(lambda x: sub_map2.get(norm_id(x)))
ic50_raw = ic50_raw.dropna(subset=['G']).copy()
ic50_raw['G'] = ic50_raw['G'].replace({'Subtype_1': 'S1', 'Subtype_2': 'S2'})

drugs9 = list(ic50_raw.columns[:9])
S1_ic = ic50_raw[ic50_raw.G == 'S1'][drugs9]
S2_ic = ic50_raw[ic50_raw.G == 'S2'][drugs9]
glob_axis = ic50_raw[drugs9].mean(axis=1)

cohen_rows, rho_rows = [], []
for d in drugs9:
    a, b = S1_ic[d].values, S2_ic[d].values
    n1, n2 = len(a), len(b)
    s1, s2 = a.std(ddof=1), b.std(ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    d_coh = (a.mean() - b.mean()) / sp
    rho, _ = spearmanr(ic50_raw[d], glob_axis)
    cohen_rows.append((d, d_coh))
    rho_rows.append((d, rho))

cohen_df = pd.DataFrame(cohen_rows, columns=['Drug', 'Cohen_d'])
rho_df = pd.DataFrame(rho_rows, columns=['Drug', 'rho']).sort_values('rho')

# Residual Cohen's d: regress out the global sensitivity axis to test
# whether the subtype effect is only a uniform IC50 shift.
resid_mat = ic50_raw[drugs9].sub(glob_axis, axis=0)
S1_resid = resid_mat[ic50_raw.G == 'S1']
S2_resid = resid_mat[ic50_raw.G == 'S2']

resid_rows = []
for d in drugs9:
    a, b = S1_resid[d].values, S2_resid[d].values
    n1, n2 = len(a), len(b)
    s1, s2 = a.std(ddof=1), b.std(ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    resid_rows.append((d, (a.mean() - b.mean()) / sp))
resid_df = pd.DataFrame(resid_rows,
                        columns=['Drug', 'Cohen_d_resid']).sort_values('Cohen_d_resid')

# Unified drug order for panels A/B/C: residual Cohen's d ascending
ABC_ORDER = resid_df['Drug'].tolist()
cohen_df = cohen_df.set_index('Drug').loc[ABC_ORDER].reset_index()


def spines_four(ax):
    for sp in ['top', 'bottom', 'left', 'right']:
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_color('#444444')
        ax.spines[sp].set_linewidth(LW)


# Mann-Whitney + BH FDR across all TIDE metrics
all_tide_metrics = ['TIDE', 'Dysfunction', 'Exclusion',
                    'CD274', 'IFNG', 'CD8', 'CTL.flag', 'CTL',
                    'MDSC', 'CAF', 'TAM M2']
pvals_raw = []
for m in all_tide_metrics:
    _, pv = mannwhitneyu(tide[tide.Group == 'S1'][m],
                         tide[tide.Group == 'S2'][m])
    pvals_raw.append(pv)
_, fdr_vals_tide, _, _ = multipletests(pvals_raw, method='fdr_bh')
fdr_map = dict(zip(all_tide_metrics, fdr_vals_tide))
pval_map = dict(zip(all_tide_metrics, pvals_raw))

# Figure: 180 x 210 mm
fig = plt.figure(figsize=(7.0866, 8.2677))
fig.patch.set_facecolor('white')

ax_a = fig.add_axes([0.08, 0.69, 0.38, 0.28])
ax_b = fig.add_axes([0.57, 0.69, 0.38, 0.28])
ax_c = fig.add_axes([0.08, 0.37, 0.38, 0.25])
ax_d = fig.add_axes([0.57, 0.37, 0.38, 0.25])
ax_e = fig.add_axes([0.08, 0.06, 0.38, 0.25])
ax_f = fig.add_axes([0.57, 0.06, 0.38, 0.25])

# Panel A
onco_off = onco_off.set_index('Drug').loc[ABC_ORDER].reset_index()
drug_names = onco_off['Drug'].tolist()
delta_display = onco_off['delta_lnIC50'].tolist()
p_vals = onco_off['p'].tolist()
fdr_vals = onco_off['FDR'].tolist()

colors_b = [C2 if v > 0 else C1 for v in delta_display]
ax_a.barh(range(len(drug_names)), delta_display,
          color=colors_b, alpha=0.85, height=0.6)
for i, (d, fdr) in enumerate(zip(delta_display, fdr_vals)):
    s = star(fdr)
    if s:
        x = d + 0.03 if d >= 0 else d - 0.03
        ax_a.text(x, i, s, va='center', ha='center',
                  fontsize=FS_STAR_A, fontweight='bold')
ax_a.set_yticks(range(len(drug_names)))
ax_a.set_yticklabels(drug_names, fontweight='bold', fontsize=FS_TICK)
ax_a.axvline(x=0, color='#333', lw=LW)
ax_a.set_xlabel(r'Difference in predicted IC$_{50}$ (log scale; S1 $-$ S2)',
                fontsize=FS_LABEL, fontweight='bold')
ax_a.legend(handles=[Patch(facecolor=C1, alpha=0.85, label='S1 more sensitive'),
                     Patch(facecolor=C2, alpha=0.85, label='S2 more sensitive')],
            fontsize=FS_LEG, framealpha=0.9)
ax_a.set_title('A', loc='left', fontsize=FS_TITLE, fontweight='bold', pad=3)
ax_a.tick_params(labelsize=FS_TICK)
spines_four(ax_a)

# Panel B
ax_b.barh(range(len(cohen_df)), cohen_df['Cohen_d'],
          color=[C1 if v < 0 else C2 for v in cohen_df['Cohen_d']],
          alpha=0.85, height=0.6)
xmax_b = cohen_df['Cohen_d'].max() + 0.05
for i, v in enumerate(cohen_df['Cohen_d']):
    ax_b.text(xmax_b, i, f'{v:+.2f}',
              va='center', ha='left', fontsize=FS_VAL, fontweight='bold',
              color=C1 if v < 0 else C2)
ax_b.set_xlim(right=xmax_b + 0.25)
ax_b.axvline(x=0, color='#333', lw=LW)
ax_b.set_yticks(range(len(cohen_df)))
ax_b.set_yticklabels(cohen_df['Drug'], fontweight='bold', fontsize=FS_TICK_BC)
ax_b.set_xlabel('Between-subtype difference\nin predicted drug sensitivity',
                fontsize=FS_LABEL, fontweight='bold')
ax_b.legend(handles=[Patch(facecolor=C1, alpha=0.85, label='S1 more sensitive'),
                     Patch(facecolor=C2, alpha=0.85, label='S2 more sensitive')],
            fontsize=FS_LEG, framealpha=0.9, loc='upper left')
ax_b.set_title('B', loc='left', fontsize=FS_TITLE, fontweight='bold', pad=3)
ax_b.tick_params(labelsize=FS_TICK_BC)
spines_four(ax_b)

# Panel C
ax_c.barh(range(len(resid_df)), resid_df['Cohen_d_resid'],
          color=[C1 if v < 0 else C2 for v in resid_df['Cohen_d_resid']],
          alpha=0.85, height=0.6)
xmax_c = resid_df['Cohen_d_resid'].max() + 0.05
for i, v in enumerate(resid_df['Cohen_d_resid']):
    ax_c.text(xmax_c, i, f'{v:+.2f}',
              va='center', ha='left', fontsize=FS_VAL, fontweight='bold',
              color=C1 if v < 0 else C2)
ax_c.set_xlim(right=xmax_c + 0.25)
ax_c.axvline(x=0, color='#333', lw=LW)
ax_c.set_yticks(range(len(resid_df)))
ax_c.set_yticklabels(resid_df['Drug'], fontweight='bold', fontsize=FS_TICK_BC)
ax_c.set_xlabel('Residual difference in predicted drug sensitivity',
                fontsize=FS_LABEL, fontweight='bold')
ax_c.legend(handles=[Patch(facecolor=C1, alpha=0.85, label='S1 more sensitive'),
                     Patch(facecolor=C2, alpha=0.85, label='S2 more sensitive')],
            fontsize=FS_LEG, framealpha=0.9, loc='upper left')
ax_c.set_title('C', loc='left', fontsize=FS_TITLE, fontweight='bold', pad=3)
ax_c.tick_params(labelsize=FS_TICK_BC)
spines_four(ax_c)

# Panel D
rates = [s1_resp / s1_n * 100, s2_resp / s2_n * 100]
bars_d = ax_d.bar([0, 1], rates, color=[C1, C2],
                  edgecolor='white', lw=LW, width=0.5)
for bar, rate, n_resp, n in zip(bars_d, rates,
                                [s1_resp, s2_resp], [s1_n, s2_n]):
    ax_d.text(bar.get_x() + bar.get_width() / 2,
              bar.get_height() + 1.5,
              f'{rate:.1f}% ({n_resp}/{n})',
              ha='center', fontsize=FS_VAL, fontweight='bold', color='#333')
ax_d.set_xticks([0, 1])
ax_d.set_xticklabels(['S1', 'S2'], fontsize=FS_TICK, fontweight='bold')
ax_d.set_xlabel('Tumor subtypes', fontsize=FS_LABEL, fontweight='bold')
ax_d.set_ylabel('Predicted response rate (%)',
                fontsize=FS_LABEL, fontweight='bold')
ax_d.set_ylim(0, max(rates) + 15)
ax_d.set_title('D', loc='left', fontsize=FS_TITLE, fontweight='bold', pad=3)
ax_d.tick_params(labelsize=FS_TICK)
spines_four(ax_d)

cont_table = [[s1_resp, s1_n - s1_resp], [s2_resp, s2_n - s2_resp]]
_, fish_p = fisher_exact(cont_table)
s_d = star(fish_p)
if s_d:
    ax_d.text(0.5, max(rates) + 5, s_d, ha='center',
              fontsize=FS_STAR_D, fontweight='bold')

# Panel E
metrics_e = ['TIDE', 'Dysfunction', 'Exclusion']
y_labels_e = ['TIDE', 'Dysfunction', 'Exclusion']

for i, m in enumerate(metrics_e):
    for j, (g, c) in enumerate([('S1', C1), ('S2', C2)]):
        vals = tide[tide.Group == g][m].values
        idx = i * 3 + j
        bp_e = ax_e.boxplot([vals], positions=[idx], patch_artist=True,
                            widths=0.55,
                            medianprops=dict(color='#333', lw=LW),
                            flierprops=dict(marker='o', markersize=2.5,
                                            alpha=0.4),
                            whiskerprops=dict(lw=LW),
                            capprops=dict(lw=LW))
        bp_e['boxes'][0].set_facecolor(c)
        bp_e['boxes'][0].set_alpha(0.75)
        jitter = np.random.normal(0, 0.04, len(vals))
        ax_e.scatter(np.ones(len(vals)) * idx + jitter, vals,
                     alpha=0.15, s=3, color=c)

ax_e.set_xticks([0.5, 3.5, 6.5])
ax_e.set_xticklabels(y_labels_e, fontsize=FS_TICK, fontweight='bold')
ax_e.set_xlabel('TIDE-derived metrics', fontsize=FS_LABEL, fontweight='bold')

for i, m in enumerate(metrics_e):
    fdr = fdr_map[m]
    s = star(fdr)
    if s:
        mid = i * 3 + 0.5
        ymax = max(tide[m].max() + 0.5,
                   tide[tide.Group == 'S1'][m].max() + 0.3)
        ax_e.text(mid, ymax, s, ha='center',
                  fontsize=FS_STAR, fontweight='bold')

e_top = max(tide[m].max() for m in metrics_e) + 1.0
ax_e.set_ylim(top=e_top)
ax_e.set_ylabel('Score', fontsize=FS_LABEL, fontweight='bold')
ax_e.set_title('E', loc='left', fontsize=FS_TITLE, fontweight='bold', pad=3)
ax_e.tick_params(labelsize=FS_TICK)
spines_four(ax_e)
ax_e.axhline(y=0, color='#999', ls='--', lw=LW, alpha=0.5)

legend_e = [Patch(facecolor=C1, alpha=0.75, label='S1'),
            Patch(facecolor=C2, alpha=0.75, label='S2')]
ax_e.legend(handles=legend_e, fontsize=FS_LEG, loc='upper right',
            ncol=2, framealpha=0.9)

# Panel F
metrics_f = ['CD274', 'IFNG', 'CD8', 'CTL.flag', 'CTL',
             'MDSC', 'CAF', 'TAM M2']
labels_f = ['PD-L1', 'IFNG', 'CD8', 'CTL flag', 'CTL',
            'MDSC', 'CAF', 'TAM M2']
n_f = len(metrics_f)

for i, m in enumerate(metrics_f):
    for j, (g, c) in enumerate([('S1', C1), ('S2', C2)]):
        vals = tide[tide.Group == g][m].values
        idx = i * 3 + j
        bp_f = ax_f.boxplot([vals], positions=[idx], patch_artist=True,
                            widths=0.55,
                            medianprops=dict(color='#333', lw=LW),
                            flierprops=dict(marker='o', markersize=2.5,
                                            alpha=0.4),
                            whiskerprops=dict(lw=LW),
                            capprops=dict(lw=LW))
        bp_f['boxes'][0].set_facecolor(c)
        bp_f['boxes'][0].set_alpha(0.75)
        jitter = np.random.normal(0, 0.04, len(vals))
        ax_f.scatter(np.ones(len(vals)) * idx + jitter, vals,
                     alpha=0.15, s=3, color=c)

ax_f.set_xticks([i * 3 + 0.5 for i in range(n_f)])
ax_f.set_xticklabels(labels_f, fontsize=FS_TICK_SM, fontweight='bold',
                     rotation=15, ha='right')

for i, m in enumerate(metrics_f):
    fdr = fdr_map[m]
    s = star(fdr)
    if s:
        mid = i * 3 + 0.5
        ymax = max(tide[m].max() + 0.5,
                   tide[tide.Group == 'S1'][m].max() + 0.3,
                   tide[tide.Group == 'S2'][m].max() + 0.3)
        ax_f.text(mid, ymax, s, ha='center',
                  fontsize=FS_STAR, fontweight='bold')

ax_f.set_title('F', loc='left', fontsize=FS_TITLE, fontweight='bold', pad=3)
ax_f.tick_params(labelsize=FS_TICK_SM)
spines_four(ax_f)
ax_f.axhline(y=0, color='#999', ls='--', lw=LW, alpha=0.5)
ax_f.set_ylabel('Normalized expression / score',
                fontsize=FS_LABEL, fontweight='bold')
ax_f.set_xlabel('Immune-related features', fontsize=FS_LABEL, fontweight='bold')

legend_f = [Patch(facecolor=C1, alpha=0.75, label='S1'),
            Patch(facecolor=C2, alpha=0.75, label='S2')]
ax_f.legend(handles=legend_f, fontsize=FS_LEG, loc='upper right',
            ncol=2, framealpha=0.9)

out_png = os.path.join(FIG_OUT, 'Figure_7_therapeutic_analysis_pubsize_v2.png')
out_svg = os.path.join(FIG_OUT, 'Figure_7_therapeutic_analysis_pubsize_v2.svg')
fig.savefig(out_png, dpi=600, bbox_inches='tight', pad_inches=0.02,
            facecolor='white', edgecolor='none')
fig.savefig(out_svg, dpi=600, bbox_inches='tight', pad_inches=0.02,
            facecolor='white', edgecolor='none',
            metadata={'Date': '2026-08-12T00:00:00'})
plt.close()

try:
    import PIL.Image
    with PIL.Image.open(out_png) as im:
        w_px, h_px = im.size
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')
    print(f'PNG size: {w_px}x{h_px} px @ 600 dpi -> '
          f'{w_px/600:.2f} x {h_px/600:.2f} inch '
          f'({w_px/600*25.4:.1f} x {h_px/600*25.4:.1f} mm)')
except Exception:
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')

print('\n=== Figure 7 stats ===')
print(f'Responder: S1={s1_resp/s1_n*100:.1f}% ({s1_resp}/{s1_n}), '
      f'S2={s2_resp/s2_n*100:.1f}% ({s2_resp}/{s2_n})')
for m in all_tide_metrics:
    s = star(fdr_map[m])
    ns_marker = ' (NS)' if s == '' else f' {s}'
    print(f'{m:<10} S1={tide[tide.Group=="S1"][m].mean():+.4f}, '
          f'S2={tide[tide.Group=="S2"][m].mean():+.4f}, '
          f'P={pval_map[m]:.4f}, FDR={fdr_map[m]:.4f}{ns_marker}')

print('\n=== Panel A oncoPredict directions ===')
for d, delta, fdr in zip(drug_names, delta_display, fdr_vals):
    sens = 'S1 (red)' if delta < 0 else 'S2 (blue)'
    print(f'{d:<12} delta(S1-S2)={delta:+.4f}  {sens:<10} '
          f'FDR={fdr:.4f} {star(fdr)}')

print("\n=== Panel B Cohen's d ===")
for _, r in cohen_df.iterrows():
    print(f'{r["Drug"]:<12} Cohen_d={r["Cohen_d"]:+.3f}')

print("\n=== Panel C residual Cohen's d ===")
for _, r in resid_df.iterrows():
    bias = 'S2' if r['Cohen_d_resid'] > 0 else 'S1'
    print(f'{r["Drug"]:<12} residual_d={r["Cohen_d_resid"]:+.3f}  '
          f'({bias}-biased)')