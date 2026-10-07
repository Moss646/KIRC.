# Fig4_plot_v31_pubsize.py
# Figure 4 (genomic), publication-size version.
#   A  mutation waterfall
#   B  14-gene CNA bar chart
#   C  16-gene CNA bar chart, using the pooled-BH v2 table
# Output: results/figures/Figure_4_genomic_v31_pubsize.{png,svg}

import sys, os, io, atexit, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

warnings.filterwarnings('ignore')


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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig4_params as P

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

FS_PANEL = 10
FS_LABEL = 8
FS_TICK = 7
FS_LEG = 7
FS_LEG_T = 7
FS_STAR = 7
SPINE_LW = 0.8

DISPLAY = {}
GROUP_ORDER = ['tumor-intrinsic', 'inflammatory PLA2']

CNV16_CSV = os.path.join(P.PROJ_ROOT, 'results', 'tables',
                         'cnv_lipid_remodeling_16genes_check_v2.csv')
OUT_PNG = os.path.join(P.OUT_DIR, 'Figure_4_genomic_v31_pubsize.png')
OUT_SVG = os.path.join(P.OUT_DIR, 'Figure_4_genomic_v31_pubsize.svg')


def draw_cna_bars(ax, df, gene_span):
    BAR_H = 0.12
    y_offsets = np.array([0.35, 0.12, -0.12, -0.35])
    n = len(df)
    y = (n - 1 - np.arange(n)) * gene_span

    ax.barh(y + y_offsets[0], df['S1_del_pct'], height=BAR_H,
            color=P.DEL_CLR, alpha=0.45)
    ax.barh(y + y_offsets[1], df['S2_del_pct'], height=BAR_H,
            color=P.DEL_CLR, alpha=0.85)
    ax.barh(y + y_offsets[2], df['S1_amp_pct'], height=BAR_H,
            color=P.AMP_CLR, alpha=0.45)
    ax.barh(y + y_offsets[3], df['S2_amp_pct'], height=BAR_H,
            color=P.AMP_CLR, alpha=0.85)

    for i, row in df.iterrows():
        s_del = P.pval_to_stars(row['p_del_adj'])
        if s_del:
            ax.text(row['S1_del_pct'] + 0.8, y[i] + y_offsets[0], s_del,
                    va='center', ha='left', fontsize=FS_STAR,
                    fontweight='bold', color=P.DEL_CLR)
            ax.text(row['S2_del_pct'] + 0.8, y[i] + y_offsets[1], s_del,
                    va='center', ha='left', fontsize=FS_STAR,
                    fontweight='bold', color=P.DEL_CLR)
        s_amp = P.pval_to_stars(row['p_amp_adj'])
        if s_amp:
            ax.text(row['S1_amp_pct'] + 0.8, y[i] + y_offsets[2], s_amp,
                    va='center', ha='left', fontsize=FS_STAR,
                    fontweight='bold', color=P.AMP_CLR)
            ax.text(row['S2_amp_pct'] + 0.8, y[i] + y_offsets[3], s_amp,
                    va='center', ha='left', fontsize=FS_STAR,
                    fontweight='bold', color=P.AMP_CLR)
    return y


def add_group_separators(ax, df, y, group_order):
    idx = 0
    n = len(df)
    for grp in group_order:
        grp_genes = [r['gene'] for _, r in df.iterrows() if r['group'] == grp]
        if grp_genes:
            end_idx = idx + len(grp_genes) - 1
            if end_idx < n - 1:
                sep_y = (y[end_idx] + y[end_idx + 1]) / 2
                ax.axhline(y=sep_y, color='#999', lw=0.6, ls='--', alpha=0.5)
            idx = end_idx + 1


def main():
    mut_stats = pd.read_csv(os.path.join(P.PROC_DIR, 'fig4_mutation_stats.csv'))
    onco = pd.read_csv(os.path.join(P.PROC_DIR, 'fig4_oncoprint.csv'))
    a_samples = pd.read_csv(os.path.join(P.PROC_DIR, 'fig4_panelA_samples.csv'))
    cnv_df = pd.read_csv(os.path.join(P.PROC_DIR, 'fig4_cnv_stats.csv'))
    cnv16 = pd.read_csv(CNV16_CSV)
    cnv16['_g'] = pd.Categorical(cnv16['group'], categories=GROUP_ORDER,
                                 ordered=True)
    cnv16 = cnv16.sort_values(['_g'], kind='stable').reset_index(drop=True)

    n_s1_mut = int((a_samples['subtype'] == 'Subtype_1').sum())
    n_s2_mut = int((a_samples['subtype'] == 'Subtype_2').sum())
    mut_pats_ordered = a_samples['patient'].tolist()
    mut_pvals_adj = dict(zip(mut_stats['gene'], mut_stats['p_bh']))

    fig = plt.figure(figsize=(7.0866, 6.8))
    fig.patch.set_facecolor(P.BG)

    # Panel A: mutation waterfall
    ax_a = fig.add_axes([0.08, 0.690, 0.86, 0.212])
    ax_a.set_facecolor('white')
    for spine in ax_a.spines.values():
        spine.set_linewidth(SPINE_LW)

    onco_map = {(r['gene'], r['patient']): r['variant_class']
                for _, r in onco.iterrows()}
    for gi, gene in enumerate(P.MUT_GENES_ORDER):
        for pi, pat in enumerate(mut_pats_ordered):
            vc = onco_map.get((gene, pat))
            if vc:
                c = P.MUT_COLORS.get(vc, P.MUT_DEFAULT_COLOR)
                rgb = tuple(int(c.lstrip('#')[i:i + 2], 16) / 255
                            for i in (0, 2, 4))
                ax_a.add_patch(plt.Rectangle((pi - 0.4, gi - 0.35), 0.8, 0.7,
                                             facecolor=rgb,
                                             edgecolor='white', lw=0.1))

    mut_labels = [f'{gene} {P.pval_to_stars(mut_pvals_adj.get(gene, 1.0))}'
                  for gene in P.MUT_GENES_ORDER]
    ax_a.set_yticks(range(len(P.MUT_GENES_ORDER)))
    ax_a.set_yticklabels(mut_labels, fontsize=FS_TICK,
                         fontstyle='italic', fontweight='bold')
    ax_a.set_ylim(len(P.MUT_GENES_ORDER) - 0.5, -0.5)
    ax_a.set_xlim(-0.5, len(mut_pats_ordered) - 0.5)
    ax_a.set_xticks([])

    s1_end = n_s1_mut
    ax_a.axvline(x=s1_end - 0.5, color='black', lw=1.0, zorder=10)

    s1_x_a = 0.08 + 0.86 * (s1_end / len(mut_pats_ordered))
    s2_x_a = 0.08 + 0.86 * ((s1_end + n_s2_mut / 2) / len(mut_pats_ordered))
    fig.text(s1_x_a, 0.673, 'S1', ha='center', fontsize=FS_TICK,
             fontweight='bold', color='black')
    fig.text(s2_x_a, 0.673, 'S2', ha='center', fontsize=FS_TICK,
             fontweight='bold', color='black')

    legend_a = [Patch(facecolor=c, label=lab, edgecolor='white', lw=0.3)
                for lab, c in P.MUT_LEGEND]
    ax_a.legend(handles=legend_a, fontsize=FS_LEG, loc='lower right',
                bbox_to_anchor=(1.0, 1.02), framealpha=0.9, ncol=5,
                title='Mutation', title_fontsize=FS_LEG_T,
                handlelength=1.0, handleheight=0.8, borderaxespad=0.3)
    ax_a.set_title('A', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)

    # Panels B and C share the same row density across different gene counts
    U_COMMON = (16 - 1) * 1.40 + 1.60
    SPAN_B = (U_COMMON - 1.60) / (len(cnv_df) - 1)
    SPAN_C = 1.40
    Y0, HGT = 0.130, 0.515
    legend_b = [Patch(facecolor=P.DEL_CLR, alpha=0.45, label='S1 deletion'),
                Patch(facecolor=P.DEL_CLR, alpha=0.85, label='S2 deletion'),
                Patch(facecolor=P.AMP_CLR, alpha=0.45, label='S1 gain'),
                Patch(facecolor=P.AMP_CLR, alpha=0.85, label='S2 gain')]

    # Panel B: 14-gene CNA bar chart
    ax_b = fig.add_axes([0.08, Y0, 0.41, HGT])
    ax_b.set_facecolor('white')
    for spine in ax_b.spines.values():
        spine.set_linewidth(SPINE_LW)

    y_b = draw_cna_bars(ax_b, cnv_df, SPAN_B)
    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels([r['gene'] for _, r in cnv_df.iterrows()],
                         fontsize=FS_TICK, fontstyle='italic',
                         fontweight='bold')
    for tick_label in ax_b.get_yticklabels():
        tick_label.set_color('black')
    add_group_separators(ax_b, cnv_df, y_b, list(P.CNV_GROUPS.keys()))

    ax_b.set_xlabel('Patients with CNA (%)', fontsize=FS_LABEL,
                    fontweight='bold')
    ax_b.set_title('B', fontsize=FS_PANEL, fontweight='bold',
                   loc='left', pad=4)
    ax_b.tick_params(axis='both', labelsize=FS_TICK)
    ax_b.legend(handles=legend_b, fontsize=FS_LEG, loc='upper right',
                framealpha=0.9, ncol=2, title='CNA',
                title_fontsize=FS_LEG_T,
                handlelength=1.0, handleheight=0.8, borderaxespad=0.3)

    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(y_b[-1] - 0.80, y_b[0] + 0.80)
    ax_b.grid(True, alpha=0.10, axis='x')

    # Panel C: 16-gene CNA bar chart
    ax_c = fig.add_axes([0.58, Y0, 0.36, HGT])
    ax_c.set_facecolor('white')
    for spine in ax_c.spines.values():
        spine.set_linewidth(SPINE_LW)

    y_c = draw_cna_bars(ax_c, cnv16, SPAN_C)
    labels16 = [DISPLAY.get(r['gene'], r['gene']) for _, r in cnv16.iterrows()]
    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels(labels16, fontsize=FS_TICK, fontstyle='italic',
                         fontweight='bold')
    for tick_label in ax_c.get_yticklabels():
        tick_label.set_color('black')
    add_group_separators(ax_c, cnv16, y_c, GROUP_ORDER)

    ax_c.set_xlabel('Patients with CNA (%)', fontsize=FS_LABEL,
                    fontweight='bold')
    ax_c.set_title('C', fontsize=FS_PANEL, fontweight='bold',
                   loc='left', pad=4)
    ax_c.tick_params(axis='both', labelsize=FS_TICK)

    ax_c.set_xlim(0, 50)
    ax_c.set_ylim(y_c[-1] - 0.80, y_c[0] + 0.80)
    ax_c.grid(True, alpha=0.10, axis='x')

    fig.savefig(OUT_PNG, dpi=600, bbox_inches='tight',
                facecolor='white', pad_inches=0.05)
    fig.savefig(OUT_SVG, format='svg', bbox_inches='tight',
                facecolor='white', pad_inches=0.05)
    plt.close()
    print(f'Saved: {OUT_PNG}')
    print(f'Saved: {OUT_SVG}')


if __name__ == '__main__':
    main()