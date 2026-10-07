#!/usr/bin/env python3
# Fig3_main_v18.py
# Figure 3: transcriptomic and CPTAC proteomic characterization.
# 3 x 3 layout:
#   A volcano | B PCA+PERMANOVA | C heatmap (top-30 lipid DEGs)
#   D GSEA    | (col2 empty)   | E GSVA
#   F key 14 lipid genes (TCGA) | G CPTAC protein | H 16 lipid-remodeling genes
# Output: results/figures/Figure_3_transcriptomic_CPTAC_v18.{png,svg}

import sys, os, io, atexit
import json
import warnings
import textwrap

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu, ttest_rel
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.utils import resample

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, Ellipse

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

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_ROOT = os.path.join(PROJ_ROOT, 'data')

RAW_DIR = os.path.join(DATA_ROOT, 'raw')
PROC_DIR = os.path.join(DATA_ROOT, 'processed')
GENE_SETS_DIR = os.path.join(DATA_ROOT, 'gene_sets')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

os.environ['MPLCONFIGDIR'] = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '.matplotlib')
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 0.8

C1 = '#C65A5A'
C2 = '#4C78A8'
COLOR_DOWN = '#1565C0'
COLOR_UP = '#E53935'

KEY_GENES_14 = ['CPT1A', 'ACOX1', 'CPT2', 'ACADSB', 'ACADM', 'ACAA2',
                'CD36', 'SLC27A2', 'FASN', 'SCD', 'PLIN2',
                'PPARG', 'PPARGC1A', 'ITPKA']


def fdr_to_stars(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


# ---- Transcriptomic data ----
print('Loading transcriptomic data...')
deg_limma = pd.read_csv(os.path.join(PROC_DIR, 'deg_limma_voom.csv')).dropna(subset=['logFC'])
expr = pd.read_csv(os.path.join(RAW_DIR, 'KIRC_expr_log2_tpm.csv'), index_col=0)
sub = pd.read_csv(os.path.join(PROC_DIR, 'subtype_assignment_balanced.csv'))
gsea_df = pd.read_csv(os.path.join(PROC_DIR, 'gsea_full_results.csv'))
gsva_df = pd.read_csv(os.path.join(PROC_DIR, 'gsva_subtype_statistics.csv'))
gsva = pd.read_csv(os.path.join(PROC_DIR, 'gsva_scores.csv'), index_col=0)

sub_map = dict(zip(sub['Patient'], sub['Subtype']))
s1_pats = [p for p, v in sub_map.items() if v == 'Subtype_1']
s2_pats = [p for p, v in sub_map.items() if v == 'Subtype_2']
s1_gsva = gsva[[c for c in gsva.columns if c in s1_pats]]
s2_gsva = gsva[[c for c in gsva.columns if c in s2_pats]]
common_samples = [c for c in expr.columns if c in s1_pats or c in s2_pats]
print(f'  Samples: {len(common_samples)} ({len(s1_pats)} S1 + {len(s2_pats)} S2)')

deg_limma_sym = pd.read_csv(os.path.join(PROC_DIR, 'deg_limma_voom_with_symbols.csv'))
deg_sig_sym = deg_limma_sym[deg_limma_sym['category'].isin(['UP_in_S1', 'DOWN_in_S1'])].dropna(subset=['gene_symbol'])
deg_sig_sym = deg_sig_sym[deg_sig_sym['gene_symbol'] != '']

deg_ensg_all = deg_limma[deg_limma['category'].isin(['UP_in_S1', 'DOWN_in_S1'])]['gene'].tolist()
_counts = pd.read_csv(os.path.join(RAW_DIR, 'KIRC_counts_for_limma.csv'), index_col=0)
_counts = _counts.reindex(columns=common_samples)
_cpm = _counts.div(_counts.sum(axis=0), axis=1) * 1e6
expr_ensg = np.log2(_cpm + 1)
deg_genes = [g for g in deg_ensg_all if g in expr_ensg.index]
print(f'  DEGs for PCA/PERMANOVA: {len(deg_genes)} genes (FDR<0.05, |log2FC|>0.5)')

_lipid_gmt = os.path.join(GENE_SETS_DIR, 'msigdb_lipid_dedup_236.gmt')
LIPID_UNIVERSE = set()
if os.path.exists(_lipid_gmt):
    with open(_lipid_gmt) as _gf:
        for _line in _gf:
            _parts = _line.rstrip('\n').split('\t')
            if len(_parts) > 2:
                LIPID_UNIVERSE.update(p.strip() for p in _parts[2:] if p.strip())
print(f'  Lipid gene universe: {len(LIPID_UNIVERSE)} genes')

# ---- CPTAC data ----
print('Loading CPTAC data...')
rnaseq_cptac = pd.read_csv(os.path.join(RAW_DIR, 'rnaseq_tumor.txt'), sep='\t', index_col=0)
prot_tumor = pd.read_csv(os.path.join(RAW_DIR, 'proteome_tumor.txt'), sep='\t', index_col=0)
prot_normal = pd.read_csv(os.path.join(RAW_DIR, 'proteome_normal.txt'), sep='\t', index_col=0)

with open(os.path.join(PROC_DIR, 'subtype_classifier.json')) as f:
    clf = json.load(f)
avail_clf_genes = [g for g in clf['genes'] if g in rnaseq_cptac.index]
avail_clf_idx = [clf['genes'].index(g) for g in avail_clf_genes]
expr_clf = rnaseq_cptac.loc[avail_clf_genes].values.T
expr_sc = (expr_clf - expr_clf.mean(axis=0)) / (expr_clf.std(axis=0, ddof=1) + 1e-8)
d1 = np.sqrt(((expr_sc - np.array(clf['centroid_s1'])[avail_clf_idx]) ** 2).sum(axis=1))
d2 = np.sqrt(((expr_sc - np.array(clf['centroid_s2'])[avail_clf_idx]) ** 2).sum(axis=1))
cptac_labels = np.where(d1 < d2, 0, 1)
cptac_samps = list(rnaseq_cptac.columns)
cptac_s1 = np.array(cptac_samps)[cptac_labels == 0]
cptac_s2 = np.array(cptac_samps)[cptac_labels == 1]
cptac_matched = sorted(set(prot_tumor.columns) & set(prot_normal.columns))
print(f'  CPTAC: {len(cptac_samps)} samples ({len(cptac_s1)} S1 + {len(cptac_s2)} S2)')

cptac_key_avail = [g for g in KEY_GENES_14 if g in prot_tumor.index]
print(f'  CPTAC key genes available: {len(cptac_key_avail)}/{len(KEY_GENES_14)}')

# ---- CPTAC S1 vs S2 (Mann-Whitney) ----
sub_data = []
for g in cptac_key_avail:
    s1 = prot_tumor.loc[g, cptac_s1].values.astype(float)
    s2 = prot_tumor.loc[g, cptac_s2].values.astype(float)
    s1 = s1[~np.isnan(s1)]
    s2 = s2[~np.isnan(s2)]
    if len(s1) < 3 or len(s2) < 3:
        continue
    _, p = mannwhitneyu(s1, s2, alternative='two-sided')
    m1, m2 = s1.mean(), s2.mean()
    se = np.sqrt(s1.var(ddof=1) / len(s1) + s2.var(ddof=1) / len(s2))
    sub_data.append({'gene': g, 'log2FC': m1 - m2,
                     'ci_low': (m1 - m2) - 1.96 * se,
                     'ci_high': (m1 - m2) + 1.96 * se,
                     'p': p, 'n_s1': len(s1), 'n_s2': len(s2)})
sub_df = pd.DataFrame(sub_data).sort_values('log2FC')
sub_df['FDR'] = multipletests(sub_df['p'].values, method='fdr_bh')[1]

tn_data = []
for g in cptac_key_avail:
    t = prot_tumor.loc[g, list(cptac_matched)].values.astype(float)
    n = prot_normal.loc[g, list(cptac_matched)].values.astype(float)
    v = ~(np.isnan(t) | np.isnan(n))
    if v.sum() < 3:
        continue
    diff = t[v] - n[v]
    m = diff.mean()
    se = diff.std() / np.sqrt(len(diff))
    _, p = ttest_rel(t[v], n[v])
    tn_data.append({'gene': g, 'log2FC': m, 'ci_low': m - 1.96 * se,
                    'ci_high': m + 1.96 * se, 'p': p, 'n': len(diff)})
tn_df = pd.DataFrame(tn_data).sort_values('log2FC')

# ---- Figure ----
print('Rendering figure...')
fig = plt.figure(figsize=(7.087, 6.18))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.30, top=0.965)
gs.update(bottom=0.07, left=0.06, right=0.96)

# Panel A: volcano (limma-voom)
ax = plt.subplot(gs[0, 0]); ax.set_facecolor('white')
n_limma_up = int((deg_limma['category'] == 'UP_in_S1').sum())
n_limma_down = int((deg_limma['category'] == 'DOWN_in_S1').sum())
n_limma_ns = int((deg_limma['category'] == 'NS').sum())
nlp = -np.log10(np.maximum(deg_limma['FDR'].values, 1e-300))
for cat, col, alpha, label in [
    ('NS', '#BDBDBD', 0.18, f'NS ({n_limma_ns})'),
    ('UP_in_S1', C1, 0.75, f'Up in S1 ({n_limma_up})'),
    ('DOWN_in_S1', C2, 0.75, f'Down in S1 ({n_limma_down})'),
]:
    m = deg_limma['category'] == cat
    ax.scatter(deg_limma.loc[m, 'logFC'], nlp[m], c=col, alpha=alpha, s=2,
               edgecolors='none', label=label)
for g_name in ['NUP93', 'MYO1E', 'COL7A1', 'TGFBI', 'ACADM', 'ACADSB', 'ECHS1', 'CPT2']:
    r = deg_limma[deg_limma['gene'] == g_name]
    if len(r):
        ax.annotate(g_name, (r['logFC'].values[0], -np.log10(r['FDR'].values[0])),
                    fontsize=5, fontweight='bold', ha='center', color='black')
ax.axvline(x=0.5, color='#333', ls='--', lw=0.6, alpha=0.4)
ax.axvline(x=-0.5, color='#333', ls='--', lw=0.6, alpha=0.4)
ax.axhline(y=-np.log10(0.05), color='#333', ls='--', lw=0.6, alpha=0.4)
ax.axvline(x=0, color='#333', ls='--', lw=0.5)
ax.set_xlabel('log2FC (S1/S2)', fontsize=8)
ax.set_ylabel('-log10(FDR)', fontsize=8)
ax.text(-0.12, 1.04, 'A', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')
ax.legend(fontsize=7, loc='upper right', framealpha=0.85)
ax.tick_params(labelsize=7)

# Panel B: PCA + PERMANOVA
ax = plt.subplot(gs[0, 1]); ax.set_facecolor('white')
X_pca_in = expr_ensg.loc[deg_genes, common_samples].T.values
X_pca_in = StandardScaler().fit_transform(X_pca_in)
pca = PCA(n_components=2, random_state=42)
X_P = pca.fit_transform(X_pca_in)
colors_pca = [C1 if sub_map.get(s) == 'Subtype_1' else C2 for s in common_samples]
ax.scatter(X_P[:, 0], X_P[:, 1], c=colors_pca, alpha=0.55, s=8, edgecolors='none')
for label, color in [('Subtype_1', C1), ('Subtype_2', C2)]:
    mask = [sub_map.get(s) == label for s in common_samples]
    if sum(mask) < 3:
        continue
    xy = X_P[mask]
    mean = xy.mean(axis=0)
    cov = np.cov(xy.T)
    evals, evecs = np.linalg.eigh(cov)
    angle = np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0]))
    w, h = 2 * np.sqrt(evals) * 1.5
    ell = Ellipse(xy=mean, width=w, height=h, angle=angle,
                  facecolor=color, alpha=0.15, edgecolor=color, lw=1.0, ls='--')
    ax.add_patch(ell)
x_p1 = pca.explained_variance_ratio_[0] * 100
x_p2 = pca.explained_variance_ratio_[1] * 100
ax.set_xlabel(f'PC1 ({x_p1:.1f}%)', fontsize=8)
ax.set_ylabel(f'PC2 ({x_p2:.1f}%)', fontsize=8)
ax.text(-0.12, 1.04, 'B', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')

s1_idx = np.array([sub_map.get(s) == 'Subtype_1' for s in common_samples])
s2_idx = ~s1_idx
D_sq = squareform(pdist(X_pca_in, 'sqeuclidean'))
N = len(common_samples)
n1 = s1_idx.sum()
n2 = s2_idx.sum()
SS_total = D_sq.sum() / N
SS_within = (D_sq[s1_idx][:, s1_idx].sum() / n1
             + D_sq[s2_idx][:, s2_idx].sum() / n2)
SS_between = SS_total - SS_within
f_obs = (SS_between / (2 - 1)) / (SS_within / (N - 2))
R2 = SS_between / SS_total
fp = np.zeros(999)
for i in range(999):
    pi = resample(np.arange(N), replace=False, n_samples=N)
    p1 = pi[:n1]
    p2 = pi[n1:]
    W = D_sq[p1][:, p1].sum() / n1 + D_sq[p2][:, p2].sum() / n2
    B = D_sq.sum() / N - W
    fp[i] = (B / (2 - 1)) / (W / (N - 2))
p_perm = (np.sum(fp >= f_obs) + 1) / 1000
ax.text(0.95, 0.05, f'PERMANOVA: R\u00b2={R2:.3f}, p={p_perm:.3f}',
        transform=ax.transAxes, fontsize=6.5, ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
ax.legend(handles=[Patch(facecolor=C1, alpha=0.6), Patch(facecolor=C2, alpha=0.6)],
          labels=['S1', 'S2'], fontsize=7, loc='upper right')
ax.tick_params(labelsize=7)

# Panel C: heatmap (top-30 lipid DEGs)
ax = plt.subplot(gs[1, 0]); ax.set_facecolor('white')
_LIPID_KW = ['FATTY_ACID', 'LIPID', 'CARNITINE', 'TRIGLYCERIDE', 'PHOSPHOLIPID',
             'CHOLESTEROL', 'PEROXISOM', 'BETA_OXIDATION', 'ACYL', 'SPHINGOLIPID',
             'GLYCEROLIPID', 'GLYCEROPHOSPHOLIPID', 'ETHER_LIPID', 'UNSATURATED',
             'STEROID', 'BILE_ACID']
_PATHWAY_LABELS = {
    'REACTOME_METABOLISM_OF_LIPIDS': 'Lipid metabolism (Reactome)',
    'KEGG_FATTY_ACID_METABOLISM': 'FA metabolism (KEGG)',
    'HALLMARK_FATTY_ACID_METABOLISM': 'FA metabolism (Hallmark)',
    'HALLMARK_BILE_ACID_METABOLISM': 'Bile acid metabolism',
    'REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA': 'PPAR\u03b1 lipid regulation',
    'REACTOME_PEROXISOMAL_LIPID_METABOLISM': 'Peroxisomal lipid metab.',
    'REACTOME_FATTY_ACID_METABOLISM': 'FA metabolism (Reactome)',
    'KEGG_PEROXISOME': 'Peroxisome (KEGG)',
    'REACTOME_SPHINGOLIPID_METABOLISM': 'Sphingolipid metabolism',
    'REACTOME_METABOLISM_OF_STEROIDS': 'Steroid metabolism',
    'REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS': 'Phospholipids/phagocytosis',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PG': 'Acyl-chain remodeling (PG)',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PI': 'Acyl-chain remodeling (PI)',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PE': 'Acyl-chain remodeling (PE)',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PC': 'Acyl-chain remodeling (PC)',
    'GOBP_ACYL_COA_METABOLIC_PROCESS': 'Acyl-CoA metabolic process',
    'GOBP_REGULATION_OF_TRIGLYCERIDE_METABOLIC_PROCESS': 'Triglyceride metab. reg.',
}
_TOP_N = 5
_lipid_mask = gsea_df['pathway'].str.contains('|'.join(_LIPID_KW), case=False)
_sig_lipid = gsea_df[_lipid_mask & (gsea_df['FDR'] < 0.05)].copy()
_left = _sig_lipid[_sig_lipid['NES'] < 0].nsmallest(_TOP_N, 'FDR')
_right = _sig_lipid[_sig_lipid['NES'] > 0].nsmallest(_TOP_N, 'FDR')
left_d = [{'label': _PATHWAY_LABELS.get(r['pathway'], r['pathway']),
           'NES': r['NES'], 'FDR': r['FDR']} for _, r in _left.iterrows()]
right_d = [{'label': _PATHWAY_LABELS.get(r['pathway'], r['pathway']),
            'NES': r['NES'], 'FDR': r['FDR']} for _, r in _right.iterrows()]
all_l = [d['label'] for d in left_d] + [d['label'] for d in right_d]
all_n = [d['NES'] for d in left_d] + [d['NES'] for d in right_d]
yp = list(range(len(all_l)))
bars = ax.barh(yp, all_n,
               color=[C2 if i < len(left_d) else C1 for i in range(len(all_l))],
               alpha=0.85, height=0.55)
panel_c_stars = {'*': 0, '**': 0, '***': 0}
fdrs_c = [d['FDR'] for d in left_d] + [d['FDR'] for d in right_d]
for i, (bar, nes, fdr) in enumerate(zip(bars, all_n, fdrs_c)):
    stars = fdr_to_stars(fdr)
    if stars:
        panel_c_stars[stars] = panel_c_stars.get(stars, 0) + 1
        ax.text(nes + (0.05 if nes >= 0 else -0.05), i, stars,
                va='center', ha='left' if nes >= 0 else 'right',
                fontsize=6.5, fontweight='bold', color='red')
print(f'Panel C: {len(all_l)} pathways, stars: {panel_c_stars}')
ax.set_yticks(yp)
ax.set_yticklabels(all_l, fontsize=5.5, fontweight='bold')
ax.yaxis.set_ticks_position('right')
ax.yaxis.set_label_position('right')
ax.axvline(x=0, color='#333', lw=0.9)
sep_y = len(left_d) - 0.5
ax.axhline(y=sep_y, color='#333', lw=1.0, alpha=0.5)
xmx = max(abs(min(all_n)), abs(max(all_n))) * 1.3
ax.set_xlim(-xmx, xmx)
ax.set_xlabel('Normalized Enrichment Score', fontsize=8)
ax.text(-0.12, 1.04, 'D', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')
ax.legend(handles=[Patch(facecolor=C1, alpha=0.75), Patch(facecolor=C2, alpha=0.75)],
          labels=['S1', 'S2'], fontsize=7, loc='upper right', framealpha=0.9)
ax.tick_params(axis='x', labelsize=7)

# Panel D: heatmap of top-30 lipid DEGs
ax = plt.subplot(gs[0, 2]); ax.set_facecolor('white')
deg_lipid = deg_limma_sym[deg_limma_sym['gene_symbol'].isin(LIPID_UNIVERSE)]
deg_lipid_sig = deg_lipid[deg_lipid['FDR'] < 0.05].sort_values('FDR')
plottable = [g for g in deg_lipid_sig['gene_symbol'] if g in expr.index]
heat_genes = plottable[:30]
rank_at_30 = 30 + sum(1 for g in deg_lipid_sig['gene_symbol'].tolist()[:30]
                      if g not in expr.index)
print(f'Panel D: {len(heat_genes)} top lipid-related DEGs '
      f'(of {len(LIPID_UNIVERSE)} lipid genes); filled from top-{rank_at_30}')
pat_sorted = sorted(common_samples, key=lambda p: (sub_map.get(p, 'Z'), p))
heat_data = expr.loc[heat_genes, pat_sorted].values
heat_z = (heat_data - heat_data.mean(axis=1, keepdims=True)) / (
    heat_data.std(axis=1, keepdims=True) + 1e-8)
im = ax.imshow(heat_z, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5,
               interpolation='nearest')
ax.set_yticks(range(len(heat_genes)))
ax.set_yticklabels(heat_genes, fontsize=4.5)
ax.set_xticks([])
s1_end = sum(1 for p in pat_sorted if sub_map.get(p) == 'Subtype_1')
ax.axvline(x=s1_end - 0.5, color='black', lw=1.2)
ax.text(s1_end / 2, len(heat_genes) + 0.5, 'S1', ha='center',
        fontsize=6.5, fontweight='bold', color=C1)
ax.text(s1_end + (len(pat_sorted) - s1_end) / 2, len(heat_genes) + 0.5,
        'S2', ha='center', fontsize=6.5, fontweight='bold', color=C2)
cbar_d = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.015)
cbar_d.set_label('Z-score', fontsize=7)
ax.text(-0.08, 1.04, 'C', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')

# Panel E: GSVA box plots
ax = plt.subplot(gs[1, 2]); ax.set_facecolor('white')
lipid_kw = ['FATTY', 'LIPID', 'CARNITINE', 'TRIGLYCERIDE', 'PHOSPHOLIPID',
            'CHOLESTEROL', 'PEROXISOM', 'BETA_OXIDATION', 'ACYL', 'SPHINGOLIPID',
            'GLYCEROLIPID', 'GLYCEROPHOSPHOLIPID', 'ETHER_LIPID', 'UNSATURATED',
            'STEROID', 'BILE_ACID', 'KETONE', 'PPARA', 'PPARG', 'SREBF', 'SCD',
            'FASN', 'ACOX', 'CPT', 'CAD', 'LIPOLYSIS', 'LIPOGENESIS']
gsva_sig = gsva_df[gsva_df['FDR'] < 0.05].copy()
gsva_sig['is_lipid'] = gsva_sig['pathway'].str.contains('|'.join(lipid_kw),
                                                        case=False, na=False)
gsva_lipid = gsva_sig[gsva_sig['is_lipid']].copy().sort_values('FDR')
lipid_down = gsva_lipid[gsva_lipid['mean_diff'] < 0].head(5)
lipid_up = gsva_lipid[gsva_lipid['mean_diff'] > 0].head(5)
gsva_sel = pd.concat([lipid_down, lipid_up], ignore_index=True)
pwy_list = gsva_sel['pathway'].tolist()

_E_SHORT = {
    'GOBP_REGULATION_OF_INTRACELLULAR_STEROID_HORMONE_RECEPTOR_SIGNALING_PATHWAY': 'Steroid recept. reg.',
    'REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA': 'PPAR\u03b1 lipid regulation',
    'GOBP_STEROID_HORMONE_RECEPTOR_SIGNALING_PATHWAY': 'Steroid recept. signaling',
    'GOBP_PEROXISOME_ORGANIZATION': 'Peroxisome organization',
    'GOBP_NEGATIVE_REGULATION_OF_INTRACELLULAR_STEROID_HORMONE_RECEPTOR_SIGNALING_PATHWAY': 'Neg. reg. steroid recept.',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PG': 'Acyl-chain remodeling (PG)',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PC': 'Acyl-chain remodeling (PC)',
    'REACTOME_ACYL_CHAIN_REMODELLING_OF_PE': 'Acyl-chain remodeling (PE)',
    'GOBP_REGULATION_OF_TRIGLYCERIDE_METABOLIC_PROCESS': 'Triglyceride metab. reg.',
    'REACTOME_LYSOSPHINGOLIPID_AND_LPA_RECEPTORS': 'Lysosphingolipid/LPA rec.',
}


def format_pwy(name):
    name = name.replace('REACTOME_', '').replace('HALLMARK_', '').replace('KEGG_', '').replace('GOBP_', '')
    name = name.replace('_', ' ')
    return textwrap.fill(name, width=22)


labels_e = [_E_SHORT.get(p, format_pwy(p)) for p in pwy_list]
fdrs_e = gsva_sel['FDR'].tolist()
labels_e = [f'{lb} {fdr_to_stars(fdr)}' for lb, fdr in zip(labels_e, fdrs_e)]
panel_e_stars = {'*': 0, '**': 0, '***': 0}
for fdr in fdrs_e:
    stars = fdr_to_stars(fdr)
    if stars:
        panel_e_stars[stars] = panel_e_stars.get(stars, 0) + 1
print(f'Panel E: {len(pwy_list)} pathways, stars: {panel_e_stars}')
y_pos = np.arange(len(pwy_list))
for idx, pwy in enumerate(pwy_list):
    if pwy not in s1_gsva.index:
        continue
    v1 = s1_gsva.loc[pwy].astype(float)
    v2 = s2_gsva.loc[pwy].astype(float)
    bp2 = ax.boxplot([v2], positions=[idx + 0.18], vert=False, widths=0.3,
                     patch_artist=True, showfliers=False,
                     whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                     medianprops=dict(color='black', lw=0.9))
    bp1 = ax.boxplot([v1], positions=[idx - 0.18], vert=False, widths=0.3,
                     patch_artist=True, showfliers=False,
                     whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                     medianprops=dict(color='black', lw=0.9))
    bp2['boxes'][0].set_facecolor(C2)
    bp2['boxes'][0].set_alpha(0.75)
    bp1['boxes'][0].set_facecolor(C1)
    bp1['boxes'][0].set_alpha(0.75)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels_e, fontsize=5.5, fontweight='bold')
ax.axvline(x=0, color='#555', ls='--', lw=0.9)
ax.set_xlabel('GSVA Enrichment Score', fontsize=8)
ax.text(-0.12, 1.04, 'E', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')
ax.legend(handles=[Patch(facecolor=C1, alpha=0.75), Patch(facecolor=C2, alpha=0.75)],
          labels=['S1', 'S2'], fontsize=7, loc='upper right', framealpha=0.9)
ax.tick_params(axis='x', labelsize=7)

# Panel F: 14 key lipid gene boxplots (TCGA)
ax = plt.subplot(gs[2, 0]); ax.set_facecolor('white')
key_genes_14 = ['CPT1A', 'ACOX1', 'CPT2', 'ACADSB', 'ACADM', 'ACAA2',
                'CD36', 'SLC27A2', 'FASN', 'SCD', 'PLIN2',
                'PPARG', 'PPARGC1A', 'ITPKA']
key_avail = [g for g in key_genes_14 if g in expr.index]
print(f'Panel F: {len(key_avail)}/{len(key_genes_14)} genes available')

limma_fdr = {}
for _, row in deg_limma_sym.iterrows():
    sym = row.get('gene_symbol', '')
    if pd.notna(sym) and sym != '' and sym in key_avail:
        limma_fdr[sym] = row['FDR']
print(f'  Limma FDR matched: {len(limma_fdr)}/{len(key_avail)} genes')

panel_f_stars = {'*': 0, '**': 0, '***': 0}
s1_expr = expr.loc[key_avail, [c for c in expr.columns if c in s1_pats]]
s2_expr = expr.loc[key_avail, [c for c in expr.columns if c in s2_pats]]
for i, gene in enumerate(key_avail):
    v1 = s1_expr.loc[gene].values.astype(float)
    v2 = s2_expr.loc[gene].values.astype(float)
    bp1 = ax.boxplot([v1], positions=[i * 2 - 0.35], widths=0.5,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='black', lw=0.9))
    bp2 = ax.boxplot([v2], positions=[i * 2 + 0.35], widths=0.5,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='black', lw=0.9))
    bp1['boxes'][0].set_facecolor(C1)
    bp1['boxes'][0].set_alpha(0.75)
    bp2['boxes'][0].set_facecolor(C2)
    bp2['boxes'][0].set_alpha(0.75)
    if gene in limma_fdr:
        stars = fdr_to_stars(limma_fdr[gene])
        if stars:
            panel_f_stars[stars] = panel_f_stars.get(stars, 0) + 1
            y_max = max(max(v1), max(v2))
            ax.text(i * 2, y_max * 1.05, stars, ha='center', va='bottom',
                    fontsize=6.5, fontweight='bold', color='red')
print(f'  Panel F stars: {panel_f_stars}')

ax.set_xticks([i * 2 for i in range(len(key_avail))])
ax.set_xticklabels(key_avail, fontsize=5.5, fontweight='bold',
                   rotation=45, ha='right')
ax.set_ylabel('log2(TPM+1)', fontsize=8)
ax.text(-0.12, 1.04, 'F', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')
ax.legend(handles=[Patch(facecolor=C1, alpha=0.75), Patch(facecolor=C2, alpha=0.75)],
          labels=['S1', 'S2'], fontsize=7, loc='upper right')
ax.tick_params(labelsize=7)
for _b in (5, 10):
    _sep_x = (_b * 2 + (_b + 1) * 2) / 2.0
    ax.axvline(x=_sep_x, color='#999999', lw=0.8, ls='--', alpha=0.6, zorder=0)

# Panel G: CPTAC protein boxplots
ax = plt.subplot(gs[2, 1]); ax.set_facecolor('white')
genes_g = [g for g in KEY_GENES_14 if g in prot_tumor.index]
fdr_g = dict(zip(sub_df['gene'], sub_df['FDR']))
panel_g_stars = {'*': 0, '**': 0, '***': 0, 'NS': 0}
for i, gene in enumerate(genes_g):
    a = prot_tumor.loc[gene, cptac_s1].values.astype(float)
    b = prot_tumor.loc[gene, cptac_s2].values.astype(float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    bp1 = ax.boxplot([a], positions=[i * 2 - 0.35], widths=0.5,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='black', lw=0.9),
                     whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    bp2 = ax.boxplot([b], positions=[i * 2 + 0.35], widths=0.5,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='black', lw=0.9),
                     whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    bp1['boxes'][0].set_facecolor(C1)
    bp1['boxes'][0].set_alpha(0.75)
    bp2['boxes'][0].set_facecolor(C2)
    bp2['boxes'][0].set_alpha(0.75)
    stars = fdr_to_stars(fdr_g.get(gene, 1.0))
    if stars:
        panel_g_stars[stars] = panel_g_stars.get(stars, 0) + 1
        y_max = max(max(a), max(b))
        ax.text(i * 2, y_max * 1.05, stars, ha='center', va='bottom',
                fontsize=6.5, fontweight='bold', color='red')
    else:
        panel_g_stars['NS'] += 1
print(f'Panel G: {len(genes_g)} genes, stars: {panel_g_stars}')
ax.set_xticks([i * 2 for i in range(len(genes_g))])
ax.set_xticklabels(genes_g, fontsize=5.5, fontweight='bold',
                   rotation=45, ha='right')
ax.set_ylabel('Protein abundance (log2)', fontsize=8)
ax.text(-0.12, 1.04, 'G', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')
ax.legend(handles=[Patch(facecolor=C1, alpha=0.75), Patch(facecolor=C2, alpha=0.75)],
          labels=['S1-like', 'S2-like'], fontsize=7, loc='upper right')
ax.tick_params(labelsize=7)

# Panel H: 16 lipid-remodeling genes (TCGA)
ax = plt.subplot(gs[2, 2]); ax.set_facecolor('white')
REMODEL_16 = ['LPCAT1', 'LPCAT2', 'LPCAT3', 'LPCAT4',
              'MBOAT1', 'MBOAT2', 'LCLAT1', 'MBOAT7',
              'PLA2G4A', 'PLA2G6', 'PNPLA8', 'PLA2G4C',
              'PLA2G2A', 'PLA2G4F', 'PLA2G1B', 'PLA2G2D']
N_INTRINSIC = 12
remodel_avail = [g for g in REMODEL_16 if g in expr.index]
print(f'Panel H: {len(remodel_avail)}/{len(REMODEL_16)} genes available')

limma_fdr_h = {}
for _, row in deg_limma_sym.iterrows():
    sym = row.get('gene_symbol', '')
    if pd.notna(sym) and sym != '' and sym in remodel_avail:
        limma_fdr_h[sym] = row['FDR']
print(f'  Limma FDR matched: {len(limma_fdr_h)}/{len(remodel_avail)} genes')

panel_h_stars = {'*': 0, '**': 0, '***': 0}
s1_expr_h = expr.loc[remodel_avail, [c for c in expr.columns if c in s1_pats]]
s2_expr_h = expr.loc[remodel_avail, [c for c in expr.columns if c in s2_pats]]
for i, gene in enumerate(remodel_avail):
    v1 = s1_expr_h.loc[gene].values.astype(float)
    v2 = s2_expr_h.loc[gene].values.astype(float)
    bp1 = ax.boxplot([v1], positions=[i * 2 - 0.35], widths=0.5,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='black', lw=0.9))
    bp2 = ax.boxplot([v2], positions=[i * 2 + 0.35], widths=0.5,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='black', lw=0.9))
    bp1['boxes'][0].set_facecolor(C1)
    bp1['boxes'][0].set_alpha(0.75)
    bp2['boxes'][0].set_facecolor(C2)
    bp2['boxes'][0].set_alpha(0.75)
    if gene in limma_fdr_h:
        stars = fdr_to_stars(limma_fdr_h[gene])
        if stars:
            panel_h_stars[stars] = panel_h_stars.get(stars, 0) + 1
            y_max = max(max(v1), max(v2))
            ax.text(i * 2, y_max * 1.05, stars, ha='center', va='bottom',
                    fontsize=6.5, fontweight='bold', color='red')
print(f'  Panel H stars: {panel_h_stars}')

for _b in (N_INTRINSIC - 1,):
    _sep_x = (_b * 2 + (_b + 1) * 2) / 2.0
    ax.axvline(x=_sep_x, color='#999999', lw=0.8, ls='--', alpha=0.6, zorder=0)

ax.set_xticks([i * 2 for i in range(len(remodel_avail))])
ax.set_xticklabels(remodel_avail, fontsize=5.5, fontweight='bold',
                   rotation=45, ha='right')
ax.set_ylabel('log2(TPM+1)', fontsize=8)
ax.text(-0.12, 1.04, 'H', transform=ax.transAxes, fontsize=10,
        fontweight='bold', va='bottom')
ax.legend(handles=[Patch(facecolor=C1, alpha=0.75), Patch(facecolor=C2, alpha=0.75)],
          labels=['S1', 'S2'], fontsize=7, loc='upper right')
ax.tick_params(labelsize=7)

# ---- Save ----
output = os.path.join(OUT_DIR, 'Figure_3_transcriptomic_CPTAC_v18.png')
svg_out = os.path.join(OUT_DIR, 'Figure_3_transcriptomic_CPTAC_v18.svg')
fig.savefig(output, dpi=600, bbox_inches='tight',
            facecolor='white', edgecolor='none')
fig.savefig(svg_out, format='svg', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print(f'\nSaved: {output}')
print(f'Saved: {svg_out}')
print('\n--- Asterisk summary ---')
print(f'Panel C (GSEA):  {len(all_l)} pathways, {panel_c_stars}')
print(f'Panel E (GSVA):  {len(pwy_list)} pathways, {panel_e_stars}')
print(f'Panel F (genes): {len(key_avail)} genes, {panel_f_stars}')
print(f'Panel G (CPTAC): {len(genes_g)} genes, {panel_g_stars}')
print('Done.')