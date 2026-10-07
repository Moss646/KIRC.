# make_purity_strata_v3.py
# Build purity strata for the stratified immune analyses (companion to
# run_immune_genes_strata_limma_v3.R).
# purity proxy = 1 - xCell ImmuneScore - StromaScore over 533 TCGA-KIRC samples;
# median split -> High (> median, n=266) / Low (<= median, n=267).
# Output: data/processed/immune_purity_strata_v3.csv

import sys, os, io, atexit
import numpy as np
import pandas as pd


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
XCL = os.path.join(DATA_ROOT, 'data', 'xcell', 'xcell_scores_raw.csv')
SUB = os.path.join(DATA_ROOT, 'data', 'processed',
                   'subtype_assignment_balanced.csv')
OUT = os.path.join(PROJ_ROOT, 'data', 'processed',
                   'immune_purity_strata_v3.csv')

xcell = pd.read_csv(XCL, index_col=0)
sub = pd.read_csv(SUB)
sub_map = dict(zip(sub['Patient'], sub['Subtype']))

patients = [c for c in xcell.columns if c in sub_map]
imm = xcell.loc['ImmuneScore', patients]
strm = xcell.loc['StromaScore', patients]
purity = 1.0 - (imm + strm)

cut = purity.median()
hi = sorted(purity[purity > cut].index)
lo = sorted(purity[purity <= cut].index)
print(f'purity proxy: median={cut:.6f} '
      f'range=[{purity.min():.6f}, {purity.max():.6f}]')
print(f'All={len(patients)} High={len(hi)} Low={len(lo)}')

s1 = {p for p in patients if sub_map[p] == 'Subtype_1'}
print(f'High purity S1={len(s1 & set(hi))}, Low purity S1={len(s1 & set(lo))}')

rows = [(p, sub_map[p], 'All') for p in patients]
rows += [(p, sub_map[p], 'High') for p in hi]
rows += [(p, sub_map[p], 'Low') for p in lo]
new = pd.DataFrame(rows, columns=['patient', 'subtype', 'stratum'])

if os.path.exists(OUT):
    old = pd.read_csv(OUT)
    if (new.shape == old.shape
            and list(new.columns) == list(old.columns)
            and new.reset_index(drop=True).equals(old.reset_index(drop=True))):
        print(f'identical to existing file ({len(old)} rows) -> rewritten in place')
        new.to_csv(OUT, index=False)
    else:
        alt = OUT.replace('.csv', '_recomputed.csv')
        n_diff = -1
        if new.shape == old.shape:
            n_diff = int((new.reset_index(drop=True)
                          != old.reset_index(drop=True)).any(axis=1).sum())
        print(f'WARNING: differs from existing file ({n_diff} rows) -> '
              f'wrote {alt}; existing file untouched')
        new.to_csv(alt, index=False)
else:
    new.to_csv(OUT, index=False)
    print(f'wrote {OUT}')