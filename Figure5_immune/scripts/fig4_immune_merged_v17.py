# fig4_immune_merged_v17.py
# Immune microenvironment of ccRCC lipid-metabolism subtypes (2x2 panels):
#   A  CIBERSORT cell fractions (MWU + BH)
#   B  ssGSEA immune pathways
#   C  MHC + cytotoxic genes (limma logFC / FDR)
#   D  Checkpoint genes (limma logFC / FDR)
# Output: results/figures/fig4_immune_merged_v17.{png,svg}
#         results/tables/fig4_immune_merged_C_D_limma_v16.csv

import sys, os, io, atexit, warnings
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
DATA_ROOT = os.environ.get('MEL_DATA_ROOT', PROJ_ROOT)

TCGA_DIR = os.path.join(DATA_ROOT, 'data', 'tcga_kirc')
CONS_DIR = os.path.join(DATA_ROOT, 'data', 'processed')
PRE_DIR = os.path.join(DATA_ROOT, 'data', 'processed')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
TAB_DIR = os.path.join(PROJ_ROOT, 'results', 'tables')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

# limma-voom table with gene symbols. Can be overridden via the environment
# variable LIMMA_CSV if the upstream Figure3 repository lives elsewhere.
LIMMA_CSV = os.environ.get(
    'LIMMA_CSV',
    os.path.join(DATA_ROOT, 'data', 'processed',
                 'deg_limma_voom_with_symbols.csv')
)

os.environ['MPLCONFIGDIR'] = os.path.join(SCRIPT_DIR, '.matplotlib')

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

BG = 'white'
C1 = '#C65A5A'
C2 = '#4C78A8'

FS_PANEL = 10
FS_LABEL = 8
FS_TICK = 7
FS_LEG = 7
FS_STAR = 7
SPINE_LW = 0.8

# ---- Load data ----
expr = pd.read_csv(os.path.join(TCGA_DIR, 'KIRC_expr_log2_tpm.csv'), index_col=0)
sub = pd.read_csv(os.path.join(CONS_DIR, 'subtype_assignment_balanced.csv'))
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
s1 = [p for p, v in sub_map.items() if v == 'Subtype_1']
s2 = [p for p, v in sub_map.items() if v == 'Subtype_2']
s1c = [c for c in expr.columns if c in s1]
s2c = [c for c in expr.columns if c in s2]
print(f'[Fig4] S1={len(s1c)} S2={len(s2c)} samples')


def mwu(v1, v2):
    return mannwhitneyu(v1, v2)


# Duplicate gene symbols in the limma table: keep the row with the larger |t|
limma_all = pd.read_csv(LIMMA_CSV)
limma_all['_abst'] = limma_all['t'].abs()
limma_all = (limma_all.sort_values('_abst', ascending=False)
             .drop_duplicates('gene_symbol', keep='first')
             .set_index('gene_symbol'))
print(f'[Fig4] limma table: {len(limma_all)} unique gene symbols')

# ---- Panel A: CIBERSORT (fraction data, MWU + BH) ----
cib_full = pd.read_csv(
    os.path.join(PRE_DIR, 'cibersort_results_official_full.csv'), index_col=0)
cib = cib_full[cib_full['P-value'] < 0.05].copy()
cell_types = [c for c in cib.columns
              if c not in ['P-value', 'Correlation', 'RMSE']]
cib_cells = cib[cell_types].copy()
cib_cells['Subtype'] = cib_cells.index.map(sub_map)
cib_cells = (cib_cells.dropna(subset=['Subtype'])
             .reset_index().rename(columns={'index': 'Patient'}))

cell_data = {}
for ct in cell_types:
    v1 = cib_cells[cib_cells['Subtype'] == 'Subtype_1'][ct].astype(float)
    v2 = cib_cells[cib_cells['Subtype'] == 'Subtype_2'][ct].astype(float)
    cell_data[ct] = {'S1': v1.mean() * 100, 'S2': v2.mean() * 100,
                     'diff': (v1.mean() - v2.mean()) * 100, 'FDR': mwu(v1, v2)[1]}
cell_df = pd.DataFrame(cell_data).T
_, cell_df['FDR'], _, _ = multipletests(cell_df['FDR'], method='fdr_bh')
cell_df = cell_df.sort_values('diff')

# ---- Panel B: ssGSEA immune pathways ----
pwy_df = pd.read_csv(os.path.join(CONS_DIR, 'ssgsea_immune_stats_R.csv'))
immune_pwys = ['ALLOGRAFT REJECTION', 'IL6 JAK STAT3 SIGNALING',
               'TNFA SIGNALING via NFKB', 'INFLAMMATORY RESPONSE',
               'COMPLEMENT', 'IL2 STAT5 SIGNALING',
               'INTERFERON GAMMA RESPONSE', 'INTERFERON ALPHA RESPONSE',
               'TGF BETA SIGNALING']
pwy_df = pwy_df[pwy_df['label'].isin(immune_pwys)].copy().sort_values('diff')

_pwy_sentence = {
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
pwy_df['label'] = pwy_df['label'].map(_pwy_sentence).fillna(pwy_df['label'])

# ---- Panel C: MHC + cytotoxic (limma logFC / FDR) ----
cp_genes = ['B2M', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'TAP1',
            'HLA-DRA', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'CD74',
            'PRF1', 'GZMB', 'GZMA', 'GNLY', 'NKG7', 'IFNG']
missing = [g for g in cp_genes if g not in limma_all.index]
if missing:
    raise SystemExit(f'genes missing from limma table: {missing}')
cdf = pd.DataFrame([
    {'gene': g, 'logFC': limma_all.loc[g, 'logFC'], 'FDR': limma_all.loc[g, 'FDR']}
    for g in cp_genes
])

# ---- Panel D: checkpoint genes (limma logFC / FDR) ----
ck = ['PDCD1', 'LAG3', 'CTLA4', 'HAVCR2', 'CD274', 'TIGIT']
ckl = ['PD-1', 'LAG-3', 'CTLA-4', 'TIM-3', 'PD-L1', 'TIGIT']
missing = [g for g in ck if g not in limma_all.index]
if missing:
    raise SystemExit(f'genes missing from limma table: {missing}')
ddf = pd.DataFrame([
    {'gene': g, 'label': l,
     'logFC': limma_all.loc[g, 'logFC'], 'FDR': limma_all.loc[g, 'FDR']}
    for g, l in zip(ck, ckl)
])

stats_out = pd.concat([
    cdf.assign(panel='C'),
    ddf.assign(panel='D'),
])[['panel', 'gene', 'label', 'logFC', 'FDR']]
stats_out.to_csv(os.path.join(TAB_DIR, 'fig4_immune_merged_C_D_limma_v16.csv'),
                 index=False)
print('[Fig4] stats CSV saved: fig4_immune_merged_C_D_limma_v16.csv')

# ---- Figure ----
fig = plt.figure(figsize=(7.0866, 5.567))
fig.patch.set_facecolor(BG)
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.40,
                      top=0.94, left=0.20, right=0.96, bottom=0.10)


def draw_stars_barh(ax, df, valcol, pad_frac, xmax):
    for i, (_, r) in enumerate(df.iterrows()):
        p = r['FDR']
        star = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        if star:
            v = r[valcol]
            xpos = v + xmax * pad_frac if v > 0 else v - xmax * pad_frac * 2
            ha = 'left' if v > 0 else 'right'
            ax.text(xpos, i, star, va='center', ha=ha,
                    fontsize=FS_STAR, fontweight='bold', color='#333')


# Panel A: CIBERSORT
ax = fig.add_subplot(gs[0, 0])
ax.set_facecolor(BG)
for spine in ax.spines.values():
    spine.set_linewidth(SPINE_LW)
ax.barh(range(len(cell_df)), cell_df['diff'],
        color=[C1 if d > 0 else C2 for d in cell_df['diff']],
        alpha=0.85, height=0.6)
ax.set_yticks(range(len(cell_df)))
ax.set_yticklabels(cell_df.index, fontsize=FS_TICK)
ax.axvline(x=0, color='#333', lw=SPINE_LW)
ax.set_xlim(-4.5, 4.5)
ax.set_xlabel('Difference in immune cell fractions (S1\u2212S2)',
              fontsize=FS_LABEL, fontweight='bold')
ax.set_title('A', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
ax.legend(handles=[Patch(facecolor=C1, alpha=0.85),
                   Patch(facecolor=C2, alpha=0.85)],
          labels=['Increased in S1', 'Increased in S2'],
          fontsize=FS_LEG, loc='upper left', framealpha=0.9)
ax.tick_params(labelsize=FS_TICK)
for i, (ct, r) in enumerate(cell_df.iterrows()):
    p = r['FDR']
    star = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
    if star:
        xpos = r['diff'] + 0.1 if r['diff'] > 0 else r['diff'] - 0.15
        ha = 'left' if r['diff'] > 0 else 'right'
        ax.text(xpos, i, star, va='center', ha=ha,
                fontsize=FS_STAR, fontweight='bold', color='#333')

# Panel C: MHC + cytotoxic
ax = fig.add_subplot(gs[0, 1])
ax.set_facecolor(BG)
for spine in ax.spines.values():
    spine.set_linewidth(SPINE_LW)
ax.barh(range(len(cdf)), cdf['logFC'],
        color=[C1 if d > 0 else C2 for d in cdf['logFC']],
        alpha=0.85, height=0.6)
ax.set_yticks(range(len(cdf)))
ax.set_yticklabels(cdf['gene'], fontsize=FS_TICK, fontstyle='italic')
ax.invert_yaxis()
ax.axvline(x=0, color='#333', lw=SPINE_LW)
xmax_c = cdf['logFC'].abs().max() * 1.18
ax.set_xlim(-xmax_c, xmax_c)
ax.set_xlabel('limma log$_2$FC (S1/S2)', fontsize=FS_LABEL, fontweight='bold')
ax.set_title('C', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
ax.tick_params(labelsize=FS_TICK)
ax.legend(handles=[Patch(facecolor=C1, alpha=0.85),
                   Patch(facecolor=C2, alpha=0.85)],
          labels=['Increased in S1', 'Increased in S2'],
          fontsize=FS_LEG, loc='upper right', framealpha=0.9)
draw_stars_barh(ax, cdf, 'logFC', 0.02, xmax_c)

# Panel B: immune pathways
ax = fig.add_subplot(gs[1, 0])
ax.set_facecolor(BG)
for spine in ax.spines.values():
    spine.set_linewidth(SPINE_LW)
ax.barh(range(len(pwy_df)), pwy_df['diff'],
        color=[C1 if d > 0 else C2 for d in pwy_df['diff']],
        alpha=0.85, height=0.6)
ax.set_yticks(range(len(pwy_df)))
ax.set_yticklabels(pwy_df['label'], fontsize=FS_TICK)
ax.axvline(x=0, color='#333', lw=SPINE_LW)
for i, (_, r) in enumerate(pwy_df.iterrows()):
    p = r['p_bh']
    star = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
    if star:
        xpos = r['diff'] + 0.002 if r['diff'] > 0 else r['diff'] - 0.004
        ha = 'left' if r['diff'] > 0 else 'right'
        ax.text(xpos, i, star, va='center', ha=ha,
                fontsize=FS_STAR, fontweight='bold', color='#333')
ax.set_xlabel('Difference in normalized ssGSEA score (S1\u2212S2)',
              fontsize=FS_LABEL, fontweight='bold')
ax.set_title('B', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
ax.tick_params(labelsize=FS_TICK)
ax.legend(handles=[Patch(facecolor=C1, alpha=0.85),
                   Patch(facecolor=C2, alpha=0.85)],
          labels=['Increased in S1', 'Increased in S2'],
          fontsize=FS_LEG, loc='upper left', framealpha=0.9)

# Panel D: checkpoints
ax = fig.add_subplot(gs[1, 1])
ax.set_facecolor(BG)
for spine in ax.spines.values():
    spine.set_linewidth(SPINE_LW)
ax.barh(range(len(ddf)), ddf['logFC'],
        color=[C1 if d > 0 else C2 for d in ddf['logFC']],
        alpha=0.85, height=0.6)
ax.set_yticks(range(len(ddf)))
ax.set_yticklabels(ddf['label'], fontsize=FS_TICK, fontstyle='italic')
ax.axvline(x=0, color='#333', lw=SPINE_LW)
xmax_d = ddf['logFC'].abs().max() * 1.18
ax.set_xlim(-xmax_d, xmax_d)
ax.set_xlabel('limma log$_2$FC (S1/S2)', fontsize=FS_LABEL, fontweight='bold')
ax.set_title('D', fontsize=FS_PANEL, fontweight='bold', loc='left', pad=4)
ax.tick_params(labelsize=FS_TICK)
ax.legend(handles=[Patch(facecolor=C1, alpha=0.85),
                   Patch(facecolor=C2, alpha=0.85)],
          labels=['Increased in S1', 'Increased in S2'],
          fontsize=FS_LEG, loc='upper left', framealpha=0.9)
draw_stars_barh(ax, ddf, 'logFC', 0.02, xmax_d)

out_png = os.path.join(OUT_DIR, 'fig4_immune_merged_v17.png')
out_svg = os.path.join(OUT_DIR, 'fig4_immune_merged_v17.svg')
fig.savefig(out_png, dpi=600, facecolor=BG, edgecolor='none')
fig.savefig(out_svg, format='svg', facecolor=BG, edgecolor='none')
plt.close()

print(f'[Fig4] Saved: {out_png}')
print(f'[Fig4] Saved: {out_svg}')
print(f'[Fig4] Panel A: CIBERSORT official (P<0.05 filtered, '
      f'{len(cib_cells)} samples)')
print(f'[Fig4] Panel B: ssGSEA ({len(pwy_df)} immune pathways)')
print(f'[Fig4] Panel C: MHC + Cytotoxic - limma logFC/FDR ({len(cdf)} genes)')
print(f'[Fig4] Panel D: Checkpoints - limma logFC/FDR ({len(ddf)} genes)')
print('\n[Fig4] Panel C limma values:')
print(cdf.to_string(index=False))
print('\n[Fig4] Panel D limma values:')
print(ddf.to_string(index=False))