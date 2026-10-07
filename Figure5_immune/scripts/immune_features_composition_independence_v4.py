# immune_features_composition_independence_v4.py
# Composition independence of immune features across purity strata.
# Panels:
#   A  CIBERSORT fractions (MWU within stratum)
#   B  ssGSEA immune pathways (MWU within stratum)
#   C  MHC + cytotoxic genes (limma within stratum)
#   D  Checkpoint genes (limma within stratum)
# Purity strata: median split of 1 - xCell ImmuneScore - StromaScore.
# Output: results/figures/Supplementary_immune_features_composition_independence_v4.{png,svg}

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
TAB_DIR = os.path.join(PROJ_ROOT, 'results', 'tables')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(TAB_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

os.environ['MPLCONFIGDIR'] = os.path.join(SCRIPT_DIR, '.matplotlib')

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

BG = 'white'
INK = '#222222'
FS_PANEL = 10
FS_LABEL = 8
FS_TICK = 7
FS_LEG = 7
FS_STAR = 7
SPINE_LW = 0.8

MARKERS = {'All': ('o', '#222222', True, 0.0),
           'High purity': ('^', '#E69F00', True, -0.30),
           'Low purity': ('s', '#0072B2', True, 0.30)}
SEP = '#BBBBBB'

CP_GENES = ['B2M', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'TAP1',
            'HLA-DRA', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'CD74',
            'PRF1', 'GZMB', 'GZMA', 'GNLY', 'NKG7', 'IFNG']
CK_GENES = ['TIGIT', 'CD274', 'HAVCR2', 'CTLA4', 'LAG3', 'PDCD1']

strat_mwu = pd.read_csv(os.path.join(
    TAB_DIR, 'immune_features_composition_stratified_mwu_v2.csv'))
limma = pd.read_csv(os.path.join(
    TAB_DIR, 'immune_genes_stratified_limma_v3.csv'))
limma = limma[limma['gene_symbol'].isin(CP_GENES + CK_GENES)].copy()
limma['stratum'] = limma['stratum'].map({'All': 'All', 'High': 'High purity',
                                         'Low': 'Low purity'})


def add_stratum_points(ax, sub, order, italic):
    ypos = {f: i for i, f in enumerate(order)}
    for _, r in sub.iterrows():
        mk, col, filled, yoff = MARKERS[r['stratum']]
        y = ypos[r['feature']] + yoff
        ax.scatter(r['diff'], y, marker=mk, s=20,
                   facecolor=col if filled else BG,
                   edgecolor=col, linewidth=0.7, zorder=3)
        fdr = r['FDR']
        star = '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else ''
        if star:
            ax.annotate(star, (r['diff'], y), textcoords='offset points',
                        xytext=(2.5, 0), va='center',
                        fontsize=FS_STAR, color=INK)
    ax.axvline(0, color='#888888', lw=0.8, zorder=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=FS_TICK,
                       fontstyle='italic' if italic else 'normal')
    ax.tick_params(axis='both', labelsize=FS_TICK, length=2.5, width=0.8)
    ax.invert_yaxis()


MM = 1 / 25.4
fig, axes = plt.subplots(2, 2, figsize=(180 * MM, 158 * MM),
                         gridspec_kw={'width_ratios': [1.30, 1.0],
                                      'wspace': 0.42, 'hspace': 0.55})

for ax in axes.ravel():
    ax.set_facecolor(BG)
    for sp in ax.spines.values():
        sp.set_linewidth(SPINE_LW)

# Panel A: CIBERSORT abundance (MWU)
A_NAMES = {
    'M0': 'Macrophages M0', 'Plasma': 'Plasma cells',
    'Tregs': 'T cells regulatory (Tregs)',
    'CD4 memory activated T': 'T cells CD4 memory activated',
    'TFH': 'T cells follicular helper', 'CD8+ T': 'Cells CD8',
    'Neutrophils': 'Neutrophils', 'NK activated': 'NK cells activated',
    'NK resting': 'NK cells resting', 'Mast activated': 'Mast cells activated',
    'T gd': 'T cells gamma delta', 'B memory': 'B cells memory',
    'CD4 naive T': 'T cells CD4 naive', 'DC activated': 'Dendritic cells activated',
    'Eosinophils': 'Eosinophils', 'B naive': 'B cells naive',
    'M1': 'Macrophages M1', 'DC resting': 'Dendritic cells resting',
    'Monocytes': 'Monocytes', 'Mast resting': 'Mast cells resting',
    'M2': 'Macrophages M2', 'CD4 memory resting T': 'T cells CD4 memory resting',
}
subA = strat_mwu[strat_mwu['feature_type'] == 'Abundance (CIBERSORT)'].copy()
subA['feature'] = subA['feature'].map(A_NAMES).fillna(subA['feature'])
orderA = subA[subA['stratum'] == 'All'].sort_values(
    'diff', ascending=False)['feature'].tolist()
add_stratum_points(axes[0, 0], subA, orderA, False)
axes[0, 0].set_xlabel('Difference in immune cell fractions (S1\u2212S2, %)',
                      fontsize=FS_LABEL)
axes[0, 0].set_title('A', fontsize=FS_PANEL, fontweight='bold',
                     loc='left', pad=2)

# Panel B: ssGSEA inflammation (MWU)
B_NAMES = {
    'ALLOGRAFT REJECTION': 'Allograft rejection',
    'IL6 JAK STAT3 SIGNALING': 'IL6 JAK STAT3 signaling',
    'TNFA SIGNALING via NFKB': 'TNF\u03b1 signaling via NF-\u03baB',
    'INFLAMMATORY RESPONSE': 'Inflammatory response',
    'COMPLEMENT': 'Complement',
    'IL2 STAT5 SIGNALING': 'IL2 STAT5 signaling',
    'INTERFERON GAMMA RESPONSE': 'Interferon gamma response',
    'INTERFERON ALPHA RESPONSE': 'Interferon \u03b1 response',
    'TGF BETA SIGNALING': 'TGF-\u03b2 signaling',
}
subB = strat_mwu[strat_mwu['feature_type'] == 'Inflammation (ssGSEA)'].copy()
subB['feature'] = subB['feature'].map(B_NAMES).fillna(subB['feature'])
orderB = subB[subB['stratum'] == 'All'].sort_values(
    'diff', ascending=False)['feature'].tolist()
add_stratum_points(axes[1, 0], subB, orderB, False)
axes[1, 0].set_xlabel('Difference in normalized ssGSEA score (S1\u2212S2)',
                      fontsize=FS_LABEL)
axes[1, 0].set_title('B', fontsize=FS_PANEL, fontweight='bold',
                     loc='left', pad=2)

# Panel C: MHC + cytotoxic (limma)
subC = limma[limma['gene_symbol'].isin(CP_GENES)].copy()
subC = subC.rename(columns={'gene_symbol': 'feature', 'logFC': 'diff'})
add_stratum_points(axes[0, 1], subC, CP_GENES, True)
axes[0, 1].axhline(5.5, color=SEP, lw=0.6, zorder=1)
axes[0, 1].axhline(10.5, color=SEP, lw=0.6, zorder=1)
axes[0, 1].set_xlabel('limma log$_2$FC (S1/S2)', fontsize=FS_LABEL)
axes[0, 1].set_title('C', fontsize=FS_PANEL, fontweight='bold',
                     loc='left', pad=2)

# Panel D: checkpoints (limma)
subD = limma[limma['gene_symbol'].isin(CK_GENES)].copy()
subD = subD.rename(columns={'gene_symbol': 'feature', 'logFC': 'diff'})
add_stratum_points(axes[1, 1], subD, CK_GENES, True)
axes[1, 1].set_xlabel('limma log$_2$FC (S1/S2)', fontsize=FS_LABEL)
axes[1, 1].set_title('D', fontsize=FS_PANEL, fontweight='bold',
                     loc='left', pad=2)

# Shared legend on panel A
h = [plt.Line2D([], [], marker='o', color='#222222', linestyle='None',
                markersize=5, markerfacecolor='#222222',
                label='All samples (n = 533)'),
     plt.Line2D([], [], marker='^', color='#E69F00', linestyle='None',
                markersize=5, markerfacecolor='#E69F00',
                label='High purity (n = 266)'),
     plt.Line2D([], [], marker='s', color='#0072B2', linestyle='None',
                markersize=5, markerfacecolor='#0072B2',
                label='Low purity (n = 267)')]
axes[0, 0].legend(handles=h, fontsize=FS_LEG, frameon=False,
                  loc='upper left', handletextpad=0.3)

axes[0, 0].set_xlim(-6.8, 6.2)
for ax in (axes[1, 0], axes[0, 1], axes[1, 1]):
    xmin, xmax = ax.get_xlim()
    pad = (xmax - xmin) * 0.14
    ax.set_xlim(xmin - pad, xmax + pad)

out_png = os.path.join(
    OUT_DIR, 'Supplementary_immune_features_composition_independence_v4.png')
out_svg = os.path.join(
    OUT_DIR, 'Supplementary_immune_features_composition_independence_v4.svg')
fig.savefig(out_png, dpi=600, bbox_inches='tight', facecolor=BG, edgecolor='none')
fig.savefig(out_svg, format='svg', bbox_inches='tight', facecolor=BG, edgecolor='none')
plt.close()
print(f'[Immune-comp v4] Saved: {out_png}')
print(f'[Immune-comp v4] Saved: {out_svg}')
print(f'[Immune-comp v4] Panel C genes n = {len(CP_GENES)}, '
      f'Panel D genes n = {len(CK_GENES)}')