# 01_pseudobulk_ncc_v5.py
# Pseudobulk aggregation of GSE159115 scRNA-seq + NCC subtype classification.
# Output: results/tables/sample_subtype_mapping_v5.csv

import sys, os, io, atexit, json
import h5py
import numpy as np
import pandas as pd
import warnings

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
SCRNA_DIR = os.path.join(DATA_ROOT, 'data', 'scrnaseq')
CONS_DIR = os.path.join(DATA_ROOT, 'data', 'processed')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'tables')
os.makedirs(OUT_DIR, exist_ok=True)

CLASSIFIER = os.path.join(CONS_DIR, 'subtype_classifier.json')
OUT_CSV = os.path.join(OUT_DIR, 'sample_subtype_mapping_v5.csv')


def aggregate_sample_pseudobulk(h5_path):
    """Sum UMI across all cells in a sample, return a Series indexed by gene."""
    sample_id = h5_path.split('_SI_')[1].split('_')[0]
    sample_key = f'SI_{sample_id}'
    print(f'  Loading {sample_key} ...', end=' ')

    with h5py.File(h5_path, 'r') as f:
        g = f['GRCh38']
        gene_names = [x.decode() for x in g['gene_names'][:]]
        data = g['data'][:]
        indices = g['indices'][:]
        indptr = g['indptr'][:]
        n_genes = len(gene_names)
        n_cells = len(indptr) - 1

    counts_per_gene = np.bincount(indices, weights=data, minlength=n_genes)
    print(f'n_cells={n_cells}, n_genes={n_genes}, total_UMI={int(counts_per_gene.sum())}')

    s = pd.Series(counts_per_gene, index=gene_names)
    if s.index.duplicated().any():
        s = s.groupby(level=0).sum()
    return sample_key, s


def main():
    h5_files = sorted([f for f in os.listdir(SCRNA_DIR)
                       if f.startswith('GSM') and f.endswith('.h5')])
    print(f'Found {len(h5_files)} h5 files')

    pseudobulk = {}
    for fn in h5_files:
        sample_key, s = aggregate_sample_pseudobulk(os.path.join(SCRNA_DIR, fn))
        pseudobulk[sample_key] = s

    df = pd.DataFrame(pseudobulk).T
    print(f'\nPseudobulk matrix: {df.shape} (samples x genes)')

    libsize = df.sum(axis=1)
    print('\nLibrary sizes:')
    print(libsize.to_string())

    cpm = df.div(libsize, axis=0) * 1e6
    log_cpm = np.log2(cpm + 1)

    with open(CLASSIFIER) as f:
        clf = json.load(f)
    genes_clf = clf['genes']
    centroid_s1 = np.array(clf['centroid_s1'])
    centroid_s2 = np.array(clf['centroid_s2'])
    scaler_mean = np.array(clf['scaler_mean'])
    scaler_std = np.array(clf['scaler_std'])

    print(f'\nNCC classifier: {len(genes_clf)} genes')
    print(f'Centroids: S1={centroid_s1.shape}, S2={centroid_s2.shape}')

    available = [g for g in genes_clf if g in log_cpm.columns]
    missing = [g for g in genes_clf if g not in log_cpm.columns]
    print(f'\nClassifier genes in pseudobulk: {len(available)}/{len(genes_clf)}')
    if missing:
        print(f'Missing genes: {missing}')

    X_68 = log_cpm.reindex(columns=genes_clf)
    X_scaled = (X_68.values - scaler_mean) / scaler_std
    X_scaled_df = pd.DataFrame(X_scaled, index=X_68.index, columns=genes_clf)

    print('\nNCC subtype classification')
    rows = []
    for sample_id in X_scaled_df.index:
        x = X_scaled_df.loc[sample_id]
        mask = ~x.isna()
        if mask.sum() < len(genes_clf) * 0.5:
            print(f'  WARNING: {sample_id} has only {mask.sum()}/{len(genes_clf)} usable genes')
        x_avail = x[mask].values
        c1_avail = centroid_s1[mask.values]
        c2_avail = centroid_s2[mask.values]
        d1 = np.sqrt(np.sum((x_avail - c1_avail) ** 2))
        d2 = np.sqrt(np.sum((x_avail - c2_avail) ** 2))
        subtype = 'S1-like' if d1 < d2 else 'S2-like'
        rows.append({'sample_id': sample_id, 'subtype': subtype,
                     'dist_s1': float(d1), 'dist_s2': float(d2),
                     'n_genes_used': int(mask.sum()),
                     'libsize': int(libsize[sample_id])})
        print(f'  {sample_id}: dist_S1={d1:.3f}, dist_S2={d2:.3f} -> {subtype}'
              f'  (used {mask.sum()}/{len(genes_clf)} genes)')

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_CSV, index=False)
    print(f'\nSaved: {OUT_CSV}')
    print('\nSubtype distribution:')
    print(out_df['subtype'].value_counts().to_string())


if __name__ == '__main__':
    main()