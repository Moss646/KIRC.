#!/usr/bin/env python3
# run_tidepy.py
# Run tidepy on the TCGA-KIRC expression matrix to produce the TIDE scores
# consumed by Fig7_main.py.
#
# Input:  data/raw/kirc_expr_matched.csv  (genes x samples, log2 TPM)
# Output: data/processed/TIDE_official_results.csv

import sys, os, io, atexit
import pandas as pd
import tidepy.pred


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
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_RAW = os.path.join(REPO_ROOT, 'data', 'raw')
DATA_PROC = os.path.join(REPO_ROOT, 'data', 'processed')

expr_path = os.path.join(DATA_RAW, 'kirc_expr_matched.csv')
print(f'Loading expression matrix: {expr_path}')
expr = pd.read_csv(expr_path, index_col=0)
print(f'  Genes: {expr.shape[0]}, Samples: {expr.shape[1]}')

print("Running tidepy.pred.TIDE(cancer='Other', force_normalize=True) ...")
result = tidepy.pred.TIDE(expr, cancer='Other', force_normalize=True)

print(f'  Done. Output shape: {result.shape}')
print(f'  Columns: {list(result.columns)}')
print(f"  Responder count: {result['Responder'].sum()} / {len(result)}")
print(f"  No benefits count: {result['No benefits'].sum()} / {len(result)}")
print(f"  TIDE range: [{result['TIDE'].min():.3f}, "
      f"{result['TIDE'].max():.3f}]")

out_path = os.path.join(DATA_PROC, 'TIDE_official_results.csv')
result.to_csv(out_path)
print(f'Saved: {out_path}')
print('Done.')