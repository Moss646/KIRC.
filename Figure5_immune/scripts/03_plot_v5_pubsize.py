# 03_plot_v5_pubsize.py
# Publication-size version of the scRNA-seq supplementary figure
# (180 x 270 mm @ 600 dpi).
# Output: results/figures/Fig_supp_scrnaseq_v5_pubsize.{png,svg}

import sys, os, io, atexit, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
TABLE_DIR = os.path.join(PROJ_ROOT, 'results', 'tables')
FIG_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

CKPT_CSV = os.path.join(TABLE_DIR, 'v5_checkpoint_subtype.csv')
MHC_CSV = os.path.join(TABLE_DIR, 'v5_mhc_subtype.csv')
FIG_PNG = os.path.join(FIG_DIR, 'Fig_supp_scrnaseq_v5_pubsize.png')
FIG_SVG = os.path.join(FIG_DIR, 'Fig_supp_scrnaseq_v5_pubsize.svg')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['axes.spines.top'] = True
plt.rcParams['axes.spines.right'] = True

FS_PANEL = 10
FS_TITLE = 8
FS_LABEL = 7
FS_TICK = 7
FS_LEG = 7
SPINE_LW = 0.8

S1_COLOR = '#B71C1C'
S2_COLOR = '#0D47A1'

CELL_TYPES = ['Tumor', 'Endo_PLVAP', 'Endo_ACKR1', 'Peri', 'vSMC', 'Macro',
              'Macro_MKI67', 'Mast', 'Tcell', 'Tcell_CD8', 'Bcell', 'Plasma']
CHECKPOINT_GENES = ['PD-1', 'CTLA-4', 'LAG-3', 'TIM-3', 'TIGIT', 'PD-L1']
MHC_GENES = ['B2M', 'HLA-A', 'HLA-B', 'HLA-C']


def draw_one_subplot(ax, df, gene, val_s1, val_s2, n1_col, n2_col,
                     ylabel, show_ylabel, show_legend, ylim_pad=1.25):
    sub = df[df['gene'] == gene].set_index('celltype').reindex(CELL_TYPES)
    s1_vals = sub[val_s1].values
    s2_vals = sub[val_s2].values

    n_ct = len(CELL_TYPES)
    x = np.arange(n_ct)
    w = 0.4
    ax.bar(x - w / 2, s1_vals, w, color=S1_COLOR,
           label='S1-like' if show_legend else None,
           edgecolor='black', linewidth=0.5)
    ax.bar(x + w / 2, s2_vals, w, color=S2_COLOR,
           label='S2-like' if show_legend else None,
           edgecolor='black', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(CELL_TYPES, rotation=45, ha='right', fontsize=FS_TICK)
    ax.tick_params(axis='both', labelsize=FS_TICK, width=SPINE_LW, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_LW)
    ax.set_title(gene, fontsize=FS_TITLE, fontweight='bold')

    ymax = max(s1_vals.max(), s2_vals.max())
    ax.set_ylim(0, ymax * ylim_pad)
    if show_ylabel:
        ax.set_ylabel(ylabel, fontsize=FS_LABEL)
    if show_legend:
        ax.legend(fontsize=FS_LEG, loc='upper left', frameon=True,
                  edgecolor='black', handlelength=1.2, handleheight=0.8)


def main():
    ckpt = pd.read_csv(CKPT_CSV)
    mhc = pd.read_csv(MHC_CSV)

    MM = 1 / 25.4
    fig = plt.figure(figsize=(180 * MM, 270 * MM))
    gs = fig.add_gridspec(5, 2, wspace=0.22, hspace=0.55,
                          left=0.10, right=0.97,
                          top=0.96, bottom=0.12)

    # Panel A: checkpoint genes (3 x 2)
    pos_a = gs[0, 0].get_position(fig)
    fig.text(pos_a.x0, pos_a.y1 + 0.010, 'A',
             fontsize=FS_PANEL, fontweight='bold', va='bottom')

    for gi, gene in enumerate(CHECKPOINT_GENES):
        row = gi // 2
        col = gi % 2
        ax = fig.add_subplot(gs[row, col])
        show_legend = (gi == 0)
        draw_one_subplot(ax, ckpt, gene,
                         val_s1='S1_pct', val_s2='S2_pct',
                         n1_col='S1_n_cells', n2_col='S2_n_cells',
                         ylabel='Cells with detectable expression (%)',
                         show_ylabel=True,
                         show_legend=show_legend,
                         ylim_pad=1.25)

    # Panel B: MHC-I genes (2 x 2)
    pos_b = gs[3, 0].get_position(fig)
    fig.text(pos_b.x0, pos_b.y1 + 0.010, 'B',
             fontsize=FS_PANEL, fontweight='bold', va='bottom')

    for gi, gene in enumerate(MHC_GENES):
        row = 3 + gi // 2
        col = gi % 2
        ax = fig.add_subplot(gs[row, col])
        show_legend = (gi == 0)
        draw_one_subplot(ax, mhc, gene,
                         val_s1='S1_mean', val_s2='S2_mean',
                         n1_col='S1_n_cells', n2_col='S2_n_cells',
                         ylabel='Mean UMI count per cell',
                         show_ylabel=True,
                         show_legend=show_legend,
                         ylim_pad=1.25)

    fig.savefig(FIG_PNG, dpi=600)
    fig.savefig(FIG_SVG)
    plt.close()
    print(f'Saved: {FIG_PNG}')
    print(f'Saved: {FIG_SVG}')


if __name__ == '__main__':
    main()