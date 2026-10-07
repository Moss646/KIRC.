#!/usr/bin/env python3
# cptac_lipid_remodeling_check_v2.py
# Protein-level S1/S2 comparison for the 16 phospholipid-remodelling genes
# in the CPTAC cohort.
# Output: results/tables/cptac_lipid_remodeling_16genes_check_v2.csv

import sys, os, io, atexit, json
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

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, 'data', 'raw')
PROC = os.path.join(BASE, 'data', 'processed')
OUT = os.path.join(BASE, 'results', 'tables',
                   'cptac_lipid_remodeling_16genes_check_v2.csv')

GENES = {
    'LPCAT1':  ('LPCAT1',  'tumor-intrinsic'),
    'LPCAT2':  ('LPCAT2',  'tumor-intrinsic'),
    'LPCAT3':  ('LPCAT3',  'tumor-intrinsic'),
    'LPCAT4':  ('LPCAT4',  'tumor-intrinsic'),
    'MBOAT1':  ('MBOAT1',  'tumor-intrinsic'),
    'MBOAT2':  ('MBOAT2',  'tumor-intrinsic'),
    'LCLAT1':  ('MBOAT5 (LCLAT1)', 'tumor-intrinsic'),
    'MBOAT7':  ('MBOAT7',  'tumor-intrinsic'),
    'PLA2G4A': ('PLA2G4A', 'tumor-intrinsic'),
    'PLA2G6':  ('PLA2G6',  'tumor-intrinsic'),
    'PNPLA8':  ('PNPLA8',  'tumor-intrinsic'),
    'PLA2G4C': ('PLA2G4C', 'tumor-intrinsic'),
    'PLA2G2A': ('PLA2G2A', 'inflammatory PLA2'),
    'PLA2G4F': ('PLA2G4F', 'inflammatory PLA2'),
    'PLA2G1B': ('PLA2G1B', 'inflammatory PLA2'),
    'PLA2G2D': ('PLA2G2D', 'inflammatory PLA2'),
}

prot = pd.read_csv(os.path.join(RAW, 'proteome_tumor.txt'),
                   sep='\t', index_col=0)
rnaseq = pd.read_csv(os.path.join(RAW, 'rnaseq_tumor.txt'),
                     sep='\t', index_col=0)

# Classify CPTAC samples with the 68-gene NCC, using cohort-specific z-scores
# (same procedure as Fig3 Panel G).
with open(os.path.join(PROC, 'subtype_classifier.json')) as f:
    clf = json.load(f)

avail = [g for g in clf['genes'] if g in rnaseq.index]
idx = [clf['genes'].index(g) for g in avail]

X = rnaseq.loc[avail].values.T
X = (X - X.mean(axis=0)) / (X.std(axis=0, ddof=1) + 1e-8)
d1 = np.sqrt(((X - np.array(clf['centroid_s1'])[idx]) ** 2).sum(axis=1))
d2 = np.sqrt(((X - np.array(clf['centroid_s2'])[idx]) ** 2).sum(axis=1))

samps = list(rnaseq.columns)
s1 = np.array(samps)[d1 < d2]
s2 = np.array(samps)[d1 >= d2]
print(f'CPTAC: {len(samps)} RNA-seq samples -> {len(s1)} S1 / {len(s2)} S2')

rows, tested = [], []
for gene_id, (label, group) in GENES.items():
    if gene_id not in prot.index:
        rows.append({'gene': label, 'group': group, 'in_matrix': 'no',
                     'n_S1': np.nan, 'n_S2': np.nan,
                     'S1_mean': np.nan, 'S2_mean': np.nan,
                     'diff_S1-S2': np.nan, 'MWU_p': np.nan,
                     'BH_FDR': np.nan, 'higher': 'NA', 'sig': 'NA'})
        continue

    v1 = prot.loc[gene_id, s1].values.astype(float)
    v2 = prot.loc[gene_id, s2].values.astype(float)
    v1 = v1[~np.isnan(v1)]
    v2 = v2[~np.isnan(v2)]

    if len(v1) < 3 or len(v2) < 3:
        rows.append({'gene': label, 'group': group, 'in_matrix': 'no(low n)',
                     'n_S1': len(v1), 'n_S2': len(v2),
                     'S1_mean': np.nan, 'S2_mean': np.nan,
                     'diff_S1-S2': np.nan, 'MWU_p': np.nan,
                     'BH_FDR': np.nan, 'higher': 'NA', 'sig': 'NA'})
        continue

    p = mannwhitneyu(v1, v2, alternative='two-sided').pvalue
    tested.append((label, group, len(v1), len(v2), v1.mean(), v2.mean(), p))

res = pd.DataFrame(tested, columns=['gene', 'group', 'n_S1', 'n_S2',
                                    'S1_mean', 'S2_mean', 'MWU_p'])
res['diff_S1-S2'] = res['S1_mean'] - res['S2_mean']
res['BH_FDR'] = multipletests(res['MWU_p'].values, method='fdr_bh')[1]
res['higher'] = np.where(res['diff_S1-S2'] > 0, 'S1', 'S2')
res['sig'] = np.where(res['BH_FDR'] < 0.05, '*', 'ns')

out = pd.concat([res, pd.DataFrame([r for r in rows if r])],
                ignore_index=True)
out['in_matrix'] = out['in_matrix'].fillna('yes')
out = out.sort_values(['group', 'BH_FDR'])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
out.to_csv(OUT, index=False)

print(out.round(4).to_string(index=False))
print(f'\nSaved: {OUT}')