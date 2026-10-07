# 02_celltype_subtype_expression_v5.py
# Per-cell-type expression of checkpoint and MHC-I genes, split by NCC subtype.
# Output: results/tables/v5_checkpoint_subtype.csv, v5_mhc_subtype.csv

import sys, os, io, atexit, gzip, traceback
import warnings
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

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
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'tables')
os.makedirs(OUT_DIR, exist_ok=True)

ANNO = os.path.join(SCRNA_DIR, 'GSE159115_ccRCC_anno.csv.gz')
MAPPING = os.path.join(OUT_DIR, 'sample_subtype_mapping_v5.csv')
OUT_CKPT = os.path.join(OUT_DIR, 'v5_checkpoint_subtype.csv')
OUT_MHC = os.path.join(OUT_DIR, 'v5_mhc_subtype.csv')

CHECKPOINT_GENES = ['PD-1', 'CTLA-4', 'LAG-3', 'TIM-3', 'TIGIT', 'PD-L1']
MHC_GENES = ['B2M', 'HLA-A', 'HLA-B', 'HLA-C']

DISPLAY = {'PDCD1': 'PD-1', 'CTLA4': 'CTLA-4', 'LAG3': 'LAG-3',
           'HAVCR2': 'TIM-3', 'TIGIT': 'TIGIT', 'CD274': 'PD-L1',
           'B2M': 'B2M', 'HLA-A': 'HLA-A', 'HLA-B': 'HLA-B', 'HLA-C': 'HLA-C'}

CELL_TYPES = ['Tumor', 'Endo_PLVAP', 'Endo_ACKR1', 'Peri', 'vSMC', 'Macro',
              'Macro_MKI67', 'Mast', 'Tcell', 'Tcell_CD8', 'Bcell', 'Plasma']


def merge_label(x):
    label_map = {
        '1:Tumor': 'Tumor',
        '2:Endo_PLVAP': 'Endo_PLVAP',
        '3:Endo_ACKR1': 'Endo_ACKR1',
        '4:Peri': 'Peri',
        '5:vSMC': 'vSMC',
        '6:Macro': 'Macro',
        '7:Macro_MKI67': 'Macro_MKI67',
        '8:Bcell': 'Bcell',
        '9:Tcell': 'Tcell',
        '10:Tcell_CD8': 'Tcell_CD8',
        '11:Tcell_CD8': 'Tcell_CD8',
        '12:Mast': 'Mast',
        '13:Plasma': 'Plasma',
    }
    return label_map.get(x, None)


def process_sample(fn, sample_key, subtype, anno, accum_ckpt, accum_mhc):
    h5_path = os.path.join(SCRNA_DIR, fn)
    with h5py.File(h5_path, 'r') as f:
        g = f['GRCh38']
        gene_names_raw = g['gene_names'][:]
        gene_names = [x.decode() if isinstance(x, bytes) else x
                      for x in gene_names_raw]
        data = g['data'][:]
        indices = g['indices'][:]
        indptr = g['indptr'][:]
        barcodes_raw = g['barcodes'][:]
        barcodes = [b.decode() if isinstance(b, bytes) else b
                    for b in barcodes_raw]
        n_genes = len(gene_names)
        n_cells = len(indptr) - 1

    print(f'    n_cells={n_cells}, n_genes={n_genes}', flush=True)

    target_idx = {}
    for g_sym in DISPLAY.keys():
        positions = [i for i, gn in enumerate(gene_names) if gn == g_sym]
        if positions:
            target_idx[g_sym] = positions[0]

    cell_of_k = np.repeat(np.arange(n_cells), np.diff(indptr))
    target_gene_set = set(target_idx.values())
    target_mask = np.isin(indices, list(target_gene_set))
    target_genes_filtered = indices[target_mask]
    target_cells_filtered = cell_of_k[target_mask]
    target_data_filtered = data[target_mask]

    coo = sparse.coo_matrix(
        (target_data_filtered, (target_cells_filtered, target_genes_filtered)),
        shape=(n_cells, n_genes)
    ).tocsr()

    anno_sample = anno[anno['sample'] == sample_key]
    cell_to_celltype = dict(zip(anno_sample['cell'], anno_sample['celltype']))

    for gene_sym, gene_disp in DISPLAY.items():
        if gene_sym not in target_idx:
            continue
        gi = target_idx[gene_sym]
        col = coo[:, gi]
        cell_umi = np.asarray(col.todense()).ravel().astype(np.int32)

        for c, bc in enumerate(barcodes):
            cell_key = f'{sample_key}_{bc}'
            celltype = cell_to_celltype.get(cell_key)
            if celltype is None or pd.isna(celltype):
                continue
            key = (subtype, celltype)
            if gene_disp in CHECKPOINT_GENES:
                accum_ckpt.setdefault(gene_disp, {}).setdefault(key, []).append(
                    int(cell_umi[c]))
            else:
                accum_mhc.setdefault(gene_disp, {}).setdefault(key, []).append(
                    int(cell_umi[c]))


def main():
    anno = pd.read_csv(gzip.open(ANNO, 'rt'))
    mapping = pd.read_csv(MAPPING)
    sub_map = dict(zip(mapping['sample_id'], mapping['subtype']))

    anno['celltype'] = anno['label'].apply(merge_label)
    anno['subtype'] = anno['sample'].map(sub_map)
    anno = anno[anno['celltype'].notna()].copy()

    print('\nCell counts per (celltype, subtype):')
    counts = anno.groupby(['subtype', 'celltype']).size().unstack(fill_value=0)
    print(counts.to_string())

    h5_files = sorted([f for f in os.listdir(SCRNA_DIR)
                       if f.startswith('GSM') and f.endswith('.h5')])

    accum_ckpt = {}
    accum_mhc = {}

    for fn in h5_files:
        sample_id = fn.split('_SI_')[1].split('_')[0]
        sample_key = f'SI_{sample_id}'
        if sample_key not in sub_map:
            continue
        subtype = sub_map[sample_key]
        print(f'  Processing {sample_key} ({subtype}) ...', flush=True)
        try:
            process_sample(fn, sample_key, subtype, anno, accum_ckpt, accum_mhc)
        except Exception:
            print(f'  FAILED on {sample_key}:')
            traceback.print_exc()
            raise

    print('\nCheckpoint (% of cells with UMI > 0):')
    ckpt_rows = []
    for ct in CELL_TYPES:
        for gene in CHECKPOINT_GENES:
            s1_cells = accum_ckpt.get(gene, {}).get(('S1-like', ct), [])
            s2_cells = accum_ckpt.get(gene, {}).get(('S2-like', ct), [])
            s1_n = len(s1_cells)
            s2_n = len(s2_cells)
            s1_pct = (np.sum(np.array(s1_cells) > 0) / s1_n * 100) if s1_n > 0 else 0
            s2_pct = (np.sum(np.array(s2_cells) > 0) / s2_n * 100) if s2_n > 0 else 0
            ckpt_rows.append({'celltype': ct, 'gene': gene,
                              'S1_pct': round(s1_pct, 4), 'S1_n_cells': s1_n,
                              'S2_pct': round(s2_pct, 4), 'S2_n_cells': s2_n})
    ckpt_df = pd.DataFrame(ckpt_rows)
    ckpt_df.to_csv(OUT_CKPT, index=False)
    print(ckpt_df.to_string(index=False))

    print('\nMHC (mean UMI per cell):')
    mhc_rows = []
    for ct in CELL_TYPES:
        for gene in MHC_GENES:
            s1_cells = accum_mhc.get(gene, {}).get(('S1-like', ct), [])
            s2_cells = accum_mhc.get(gene, {}).get(('S2-like', ct), [])
            s1_n = len(s1_cells)
            s2_n = len(s2_cells)
            s1_mean = float(np.mean(s1_cells)) if s1_n > 0 else 0
            s2_mean = float(np.mean(s2_cells)) if s2_n > 0 else 0
            mhc_rows.append({'celltype': ct, 'gene': gene,
                             'S1_mean': round(s1_mean, 4), 'S1_n_cells': s1_n,
                             'S2_mean': round(s2_mean, 4), 'S2_n_cells': s2_n})
    mhc_df = pd.DataFrame(mhc_rows)
    mhc_df.to_csv(OUT_MHC, index=False)
    print(mhc_df.to_string(index=False))


if __name__ == '__main__':
    main()