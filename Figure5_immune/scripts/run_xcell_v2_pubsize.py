# run_xcell_v2_pubsize.py
# Publication-size variant of run_xcell_v2.py.
# Consumes pre-computed xCell scores (R xCell v1.1.0, spillover-corrected)
# and renders a 2x2 panel at 180 mm @ 600 dpi.
# Output: results/figures/Supplementary_xCell_analysis_v17_pubsize.{png,svg}

import sys, os, io, atexit, warnings
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_ROOT = os.environ.get('MEL_DATA_ROOT', PROJ_ROOT)
TCGA_DIR = os.path.join(DATA_ROOT, 'data', 'tcga_kirc')
CONS_DIR = os.path.join(DATA_ROOT, 'data', 'processed')
XCELL_DIR = os.path.join(DATA_ROOT, 'data', 'xcell')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8

FS_PANEL = 10
FS_LABEL = 8
FS_TICK = 7
FS_LEG = 7
FS_LEG_T = 7
FS_STAR = 7
SPINE_LW = 0.8

C_S1 = '#E74C3C'
C_S2 = '#3498DB'


def _hold():
    if os.environ.get('NOPAUSE'):
        return
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.getch()
    except Exception:
        pass


def main():
    print('[xCell] Loading expression + subtype + xCell scores...')
    EXPR_TPM = pd.read_csv(os.path.join(TCGA_DIR, 'KIRC_expr_tpm.csv'),
                           index_col=0)
    print(f'  Expression matrix: {EXPR_TPM.shape[0]} genes x '
          f'{EXPR_TPM.shape[1]} samples')

    sub = pd.read_csv(os.path.join(CONS_DIR, 'subtype_assignment_balanced.csv'))
    sub_map = dict(zip(sub['Patient'], sub['Subtype']))
    s1_samples = [p for p, v in sub_map.items() if v == 'Subtype_1']
    s2_samples = [p for p, v in sub_map.items() if v == 'Subtype_2']
    print(f'  S1: {len(s1_samples)}, S2: {len(s2_samples)}')

    xcell_scores = pd.read_csv(os.path.join(XCELL_DIR, 'xcell_scores_raw.csv'),
                               index_col=0)
    print(f'  xCell raw scores: {xcell_scores.shape[0]} cell types x '
          f'{xcell_scores.shape[1]} samples')

    xcell_z = pd.read_csv(
        os.path.join(XCELL_DIR, 'xcell_scores_zscore.csv'), index_col=0).T
    print(f'  xCell z-scores: {xcell_z.shape} (samples x cell types)')

    diff_df = pd.read_csv(os.path.join(XCELL_DIR, 'xcell_subtype_diff.csv'))
    print(f'  Diff stats: {len(diff_df)} cell types, '
          f'{(diff_df["fdr"] < 0.05).sum()} significant (FDR<0.05)')

    common = [s for s in xcell_z.index if s in sub_map]
    xcell_filt = xcell_z.loc[common].astype(float)
    s1_idx = [s for s in common if sub_map[s] == 'Subtype_1']
    s2_idx = [s for s in common if sub_map[s] == 'Subtype_2']

    print('[xCell] Generating publication-size figure...')

    plot_cts = diff_df[diff_df['fdr'] < 0.05].head(30)['cell_type'].tolist()
    if len(plot_cts) < 15:
        plot_cts = diff_df.head(20)['cell_type'].tolist()

    fig = plt.figure(figsize=(7.0866, 7.0866))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.0, 1.05])
    axes = np.array([[fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])],
                     [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]])
    plt.subplots_adjust(left=0.13, right=0.93, bottom=0.11, top=0.93,
                        wspace=0.30, hspace=0.28)

    # Panel A: differential xCell scores (S1 minus S2)
    ax = axes[0, 0]
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
    plot_data = diff_df[diff_df['cell_type'].isin(plot_cts)].copy()
    plot_data['diff_s1s2'] = -plot_data['mean_zscore_diff']
    plot_data = plot_data.sort_values('diff_s1s2', ascending=True)
    colors = [C_S1 if v > 0 else C_S2 for v in plot_data['diff_s1s2']]
    ax.barh(range(len(plot_data)), plot_data['diff_s1s2'], color=colors,
            edgecolor='white', linewidth=0.3, height=0.7)
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data['cell_type'], fontsize=FS_TICK)

    _xmin, _xmax = plot_data['diff_s1s2'].min(), plot_data['diff_s1s2'].max()
    _xpad = (_xmax - _xmin) * 0.04
    for _i, (_idx, _row) in enumerate(plot_data.iterrows()):
        _fdr = _row['fdr']
        _stars = ('***' if _fdr < 0.001 else '**' if _fdr < 0.01
                  else '*' if _fdr < 0.05 else 'ns')
        _tip = _row['diff_s1s2']
        _sx = _tip + _xpad if _tip >= 0 else _tip - _xpad
        _ha = 'left' if _tip >= 0 else 'right'
        ax.text(_sx, _i, _stars, va='center', ha=_ha, fontsize=FS_STAR,
                fontweight='bold', color='black')
    ax.set_xlim(_xmin - _xpad * 2, _xmax + _xpad * 2.5)

    ax.set_xlabel('Difference in xCell score (S1 \u2212 S2)', fontsize=FS_LABEL)
    ax.set_title('A', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
    ax.axvline(0, color='black', linewidth=0.6)
    ax.tick_params(axis='both', labelsize=FS_TICK)
    legend_a = [Patch(facecolor=C_S1, alpha=0.85, label='Enriched in S1'),
                Patch(facecolor=C_S2, alpha=0.85, label='Enriched in S2')]
    ax.legend(handles=legend_a, fontsize=FS_LEG, loc='upper left',
              framealpha=0.9, handlelength=1.0, handleheight=0.8,
              borderaxespad=0.3)

    # Panel B: key immune cell boxplots (activation / APC / suppression)
    ax = axes[0, 1]
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
    panel_b_display = [
        ('Th1 cells', 'Th1 cells'),
        ('NKT', 'NKT cells'),
        ('CD8+ T-cells', 'CD8+ T cells'),
        ('B-cells', 'B cells'),
        ('DC', 'DC'),
        ('Tregs', 'Tregs'),
        ('Macrophages', 'Macrophages'),
    ]
    panel_b_plot = list(reversed(panel_b_display))
    b_cts = [d for d, _lab in panel_b_plot]
    b_labels = [lab for _d, lab in panel_b_plot]

    bp_data, bp_positions, bp_colors = [], [], []
    for i, ct in enumerate(b_cts):
        s1_vals = xcell_filt.loc[s1_idx, ct].dropna().values
        s2_vals = xcell_filt.loc[s2_idx, ct].dropna().values
        bp_data.append(s1_vals)
        bp_positions.append(i - 0.15)
        bp_colors.append(C_S1)
        bp_data.append(s2_vals)
        bp_positions.append(i + 0.15)
        bp_colors.append(C_S2)

    bp = ax.boxplot(bp_data, positions=bp_positions, vert=False,
                    patch_artist=True, widths=0.22, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.0),
                    whiskerprops=dict(linewidth=0.5),
                    capprops=dict(linewidth=0.5),
                    boxprops=dict(linewidth=0.5))
    for patch, color in zip(bp['boxes'], bp_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor('white')

    for ysep in [1.5, 3.5]:
        ax.axhline(ysep, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    def _upper_whisker(v):
        q1, q3 = np.percentile(v, [25, 75])
        return q3 + 1.5 * (q3 - q1)

    star_x = max(_upper_whisker(d) for d in bp_data if len(d) > 0) + 0.5
    for i, ct in enumerate(b_cts):
        fdr = diff_df.loc[diff_df['cell_type'] == ct, 'fdr']
        fdr = fdr.values[0] if len(fdr) else 1.0
        stars = ('***' if fdr < 0.001 else '**' if fdr < 0.01
                 else '*' if fdr < 0.05 else 'ns')
        ax.text(star_x, i, stars, va='center', ha='left', fontsize=FS_STAR,
                fontweight='bold', color='black')

    ax.set_yticks(range(len(b_cts)))
    ax.set_yticklabels(b_labels, fontsize=FS_TICK)
    ax.set_xlabel('xCell score (row-wise Z-score)', fontsize=FS_LABEL)
    ax.set_title('B', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
    ax.set_ylim(-0.8, len(b_cts) - 0.2)
    ax.set_xlim(-3, 4)
    ax.set_xticks(np.arange(-2, 5, 2))
    ax.tick_params(axis='both', labelsize=FS_TICK)
    legend_elements = [Patch(facecolor=C_S1, alpha=0.75, label='S1'),
                       Patch(facecolor=C_S2, alpha=0.75, label='S2')]
    ax.legend(handles=legend_elements, fontsize=FS_LEG, loc='upper right',
              ncol=2, handlelength=1.0, handleheight=0.8, borderaxespad=0.3)

    # Panel C: xCell composite ImmuneScore vs StromaScore
    ax = axes[1, 0]
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
    comp = xcell_scores.loc[['StromaScore', 'ImmuneScore']]
    s1c = [c for c in s1_samples if c in comp.columns]
    s2c = [c for c in s2_samples if c in comp.columns]
    ax.scatter(comp.loc['StromaScore', s1c], comp.loc['ImmuneScore', s1c],
               c=C_S1, label='S1', alpha=0.5, s=12, marker='o',
               edgecolors='none')
    ax.scatter(comp.loc['StromaScore', s2c], comp.loc['ImmuneScore', s2c],
               c=C_S2, label='S2', alpha=0.5, s=12, marker='^',
               edgecolors='none')
    ax.set_xlabel('xCell StromaScore', fontsize=FS_LABEL)
    ax.set_ylabel('xCell ImmuneScore', fontsize=FS_LABEL)
    ax.set_title('C', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
    ax.legend(fontsize=FS_LEG, loc='upper right', markerscale=1.5,
              handlelength=1.0, borderaxespad=0.3)
    ax.tick_params(axis='both', labelsize=FS_TICK)
    immune_score_s1 = comp.loc['ImmuneScore', s1c]
    immune_score_s2 = comp.loc['ImmuneScore', s2c]

    # Panel D: heatmap of top-variance cell types
    ax = axes[1, 1]
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
    heatmap_data = pd.DataFrame({
        'S1': xcell_filt.loc[s1_idx].mean(),
        'S2': xcell_filt.loc[s2_idx].mean(),
    }).T
    variances = xcell_filt.var()
    top_var = variances.nlargest(25).index
    heatmap_data = heatmap_data[top_var]

    im = ax.imshow(heatmap_data.values, aspect='auto',
                   cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xticks(range(len(top_var)))
    ax.set_xticklabels(top_var, fontsize=FS_TICK, rotation=90, ha='center')
    ax.set_yticks(range(2))
    ax.set_yticklabels(['S1', 'S2'], fontsize=FS_TICK)
    ax.set_title('D', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
    ax.tick_params(axis='both', labelsize=FS_TICK)
    cbar = plt.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label('Mean Z-score', fontsize=FS_LABEL)
    cbar.ax.tick_params(labelsize=FS_TICK)
    cbar.outline.set_linewidth(SPINE_LW)

    out_png = os.path.join(
        OUT_DIR, 'Supplementary_xCell_analysis_v17_pubsize.png')
    out_svg = os.path.join(
        OUT_DIR, 'Supplementary_xCell_analysis_v17_pubsize.svg')
    fig.savefig(out_png, dpi=600, facecolor='white', pad_inches=0.05)
    fig.savefig(out_svg, format='svg', facecolor='white', pad_inches=0.05)
    plt.close()
    print(f'[xCell] Saved: {out_png}')
    print(f'[xCell] Saved: {out_svg}')
    print(f'[xCell] S1 mean immune score: {immune_score_s1.mean():.3f}')
    print(f'[xCell] S2 mean immune score: {immune_score_s2.mean():.3f}')
    print()
    print('Top 10 differential cell types (S2 vs S1):')
    for _, row in diff_df.head(10).iterrows():
        stars = ('***' if row['fdr'] < 0.001 else '**' if row['fdr'] < 0.01
                 else '*' if row['fdr'] < 0.05 else 'ns')
        direction = 'S2-up' if row['mean_zscore_diff'] > 0 else 'S2-down'
        print(f'  {row["cell_type"]:<30s} '
              f'mean_zscore_diff={row["mean_zscore_diff"]:+7.3f}  '
              f'FDR={row["fdr"]:.2e} {stars} {direction}')


if __name__ == '__main__':
    atexit.register(_hold)
    main()