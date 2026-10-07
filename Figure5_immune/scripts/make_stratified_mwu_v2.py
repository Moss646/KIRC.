# make_stratified_mwu_v2.py
# Stratified (purity) S1-vs-S2 immune feature table for Figure 5 panels A/B.
# Feature families: CIBERSORT abundance, ssGSEA inflammation,
# cytotoxic genes, checkpoint genes.
# Output: results/tables/immune_features_composition_stratified_mwu_v2.csv

import sys, os, io, atexit
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


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
D = os.path.join(DATA_ROOT, 'data')

CS = os.path.join(D, 'processed', 'cibersort_results_official_full.csv')
SS = os.path.join(D, 'processed', 'ssgsea_immune_scores_R.csv')
SUB = os.path.join(D, 'processed', 'subtype_assignment_balanced.csv')
STRATA = os.path.join(D, 'processed', 'immune_purity_strata_v3.csv')
XCL = os.path.join(D, 'xcell', 'xcell_scores_raw.csv')
EXPR = os.path.join(D, 'tcga_kirc', 'KIRC_expr_log2_tpm.csv')
OUT = os.path.join(PROJ_ROOT, 'results', 'tables',
                   'immune_features_composition_stratified_mwu_v2.csv')

CIBERSORT_RENAME = {
    'B cells naive': 'B naive', 'B cells memory': 'B memory',
    'Plasma cells': 'Plasma',
    'T cells CD8': 'CD8+ T', 'T cells CD4 naive': 'CD4 naive T',
    'T cells CD4 memory resting': 'CD4 memory resting T',
    'T cells CD4 memory activated': 'CD4 memory activated T',
    'T cells follicular helper': 'TFH',
    'T cells regulatory (Tregs)': 'Tregs',
    'T cells gamma delta': 'T gd', 'NK cells resting': 'NK resting',
    'NK cells activated': 'NK activated', 'Macrophages M0': 'M0',
    'Macrophages M1': 'M1', 'Macrophages M2': 'M2',
    'Dendritic cells resting': 'DC resting',
    'Dendritic cells activated': 'DC activated',
    'Mast cells resting': 'Mast resting',
    'Mast cells activated': 'Mast activated'}
SSGSEA_SETS = [
    ('ALLOGRAFT REJECTION', 'HALLMARK_ALLOGRAFT_REJECTION'),
    ('IL6 JAK STAT3 SIGNALING', 'HALLMARK_IL6_JAK_STAT3_SIGNALING'),
    ('TNFA SIGNALING via NFKB', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB'),
    ('INFLAMMATORY RESPONSE', 'HALLMARK_INFLAMMATORY_RESPONSE'),
    ('COMPLEMENT', 'HALLMARK_COMPLEMENT'),
    ('IL2 STAT5 SIGNALING', 'HALLMARK_IL2_STAT5_SIGNALING'),
    ('INTERFERON GAMMA RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
    ('INTERFERON ALPHA RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE'),
    ('TGF BETA SIGNALING', 'HALLMARK_TGF_BETA_SIGNALING')]
CYTOTOXIC = ['PRF1', 'GZMB', 'GZMA', 'GNLY', 'NKG7', 'IFNG']
CHECKPOINT = ['PDCD1', 'LAG3', 'CTLA4', 'HAVCR2', 'CD274', 'TIGIT']

sub = pd.read_csv(SUB)
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
s1 = {p for p in sub_map if sub_map[p] == 'Subtype_1'}
s2 = {p for p in sub_map if sub_map[p] == 'Subtype_2'}

strata = pd.read_csv(STRATA)
stratum_of = dict(zip(strata['patient'], strata['stratum']))


def stratum_groups(patients):
    hi = sorted(p for p in patients if stratum_of.get(p) == 'High')
    lo = sorted(p for p in patients if stratum_of.get(p) == 'Low')
    return [('All', sorted(patients)), ('High purity', hi), ('Low purity', lo)]


def mwu_rows(ftype, feature, values, patients):
    rows = []
    for lab, grp in stratum_groups(patients):
        a = values.reindex(sorted(set(grp) & s1)).dropna()
        b = values.reindex(sorted(set(grp) & s2)).dropna()
        rows.append((ftype, feature, lab, a.mean(), b.mean(),
                     a.mean() - b.mean(),
                     mannwhitneyu(a, b)[1], len(a), len(b)))
    return rows


rows = []

cs = pd.read_csv(CS)
cs = cs[cs['P-value'] < 0.05].copy()
cs_pat = cs['Unnamed: 0'].tolist()
for col in [c for c in cs.columns
            if c in list(CIBERSORT_RENAME) + ['Monocytes', 'Eosinophils',
                                              'Neutrophils']]:
    vals = cs[col] * 100.0
    vals.index = cs_pat
    rows += mwu_rows('Abundance (CIBERSORT)',
                     CIBERSORT_RENAME.get(col, col), vals, cs_pat)

ss = pd.read_csv(SS, index_col=0)
for label, set_name in SSGSEA_SETS:
    rows += mwu_rows('Inflammation (ssGSEA)', label,
                     ss.loc[set_name], list(ss.columns))

expr = pd.read_csv(EXPR, index_col=0)
for gtype, genes in [('Cytotoxic genes', CYTOTOXIC),
                     ('Checkpoint genes', CHECKPOINT)]:
    for g in genes:
        rows += mwu_rows(gtype, g, expr.loc[g], list(expr.columns))

new = pd.DataFrame(rows, columns=['feature_type', 'feature', 'stratum',
                                  'S1_mean', 'S2_mean', 'diff', 'p',
                                  'n1', 'n2'])
for (ftype, st), m in new.groupby(['feature_type', 'stratum']).groups.items():
    idx = list(m)
    new.loc[idx, 'FDR'] = multipletests(new.loc[idx, 'p'],
                                        method='fdr_bh')[1]

print(f'rows={len(new)}  CIBERSORT subset n={len(cs_pat)} '
      f'(High={sum(stratum_of.get(p) == "High" for p in cs_pat)}, '
      f'Low={sum(stratum_of.get(p) == "Low" for p in cs_pat)})')

if os.path.exists(OUT):
    old = pd.read_csv(OUT)
    if (new.shape == old.shape
            and list(new.columns) == list(old.columns)
            and np.allclose(new.select_dtypes('number'),
                            old.select_dtypes('number'), atol=1e-12)
            and new.select_dtypes('object').equals(old.select_dtypes('object'))):
        print(f'identical to existing file ({len(old)} rows) -> rewritten in place')
        new.to_csv(OUT, index=False)
    else:
        alt = OUT.replace('.csv', '_recomputed.csv')
        print(f'WARNING: differs from existing file -> wrote {alt}; '
              f'existing file untouched')
        new.to_csv(alt, index=False)
else:
    new.to_csv(OUT, index=False)
    print(f'wrote {OUT}')