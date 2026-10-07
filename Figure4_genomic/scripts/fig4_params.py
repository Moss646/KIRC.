# fig4_params.py
# Shared parameters and paths for the Figure 4 (genomic) pipeline.

import os

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_ROOT = os.path.join(PROJ_ROOT, 'data')
RAW_DIR = os.path.join(DATA_ROOT, 'raw')
PROC_DIR = os.path.join(DATA_ROOT, 'processed')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

# Panel A driver genes (top to bottom)
MUT_GENES_ORDER = ['VHL', 'PBRM1', 'BAP1', 'SETD2', 'KDM5C', 'MTOR', 'TP53']

# Nonsynonymous variant classes retained for analysis
NONSYN = ['Missense_Mutation', 'Nonsense_Mutation',
          'Frame_Shift_Del', 'Frame_Shift_Ins',
          'Splice_Site', 'In_Frame_Del', 'In_Frame_Ins']

# Fill color per variant class (Panel A oncoprint)
MUT_COLORS = {
    'Missense_Mutation': '#2E86AB',
    'Nonsense_Mutation': '#000000',
    'Frame_Shift_Del': '#F18F01',
    'Frame_Shift_Ins': '#C73E1D',
    'Splice_Site': '#8A4FFF',
    'In_Frame_Del': '#3B1F2B',
    'In_Frame_Ins': '#FFBF46',
}
MUT_DEFAULT_COLOR = '#7FB069'

MUT_LEGEND = [
    ('Missense', '#2E86AB'),
    ('Nonsense', '#000000'),
    ('Frameshift', '#F18F01'),
    ('Splice site', '#8A4FFF'),
    ('In-frame deletion/insertion', '#3B1F2B'),
]

# Panel B lipid-metabolism genes grouped by function
CNV_GROUPS = {
    'FAO': ['CPT1A', 'ACOX1', 'CPT2', 'ACADSB', 'ACADM', 'ACAA2'],
    'FA uptake/synthesis/storage': ['CD36', 'SLC27A2', 'FASN', 'SCD', 'PLIN2'],
    'Lipid regulators': ['PPARG', 'PPARGC1A', 'ITPKA'],
}
GROUP_COLORS = {
    'FAO': '#2E86AB',
    'FA uptake/synthesis/storage': '#F18F01',
    'Lipid regulators': '#8A4FFF',
}

DEL_CLR = '#B71C1C'
AMP_CLR = '#0D47A1'

BG = '#FAFAFA'


def bh_correct(pvals):
    import numpy as np
    arr = np.asarray(pvals, dtype=float)
    n = len(arr)
    if n == 0:
        return arr
    order = np.argsort(arr)
    sorted_p = arr[order]
    adj = sorted_p * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out = np.zeros(n)
    out[order] = adj
    return out


def pval_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''