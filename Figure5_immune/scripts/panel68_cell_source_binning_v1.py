# panel68_cell_source_binning_v1.py
# Bin the 68 NCC panel genes by cell source in GSE159115 scRNA-seq.
# For each gene x compartment (Tumor / Myeloid / Lymphoid / Stromal) report
# mean CP10K, % cells expressing, UMI share, and dominant compartment.
# Output: results/tables/panel68_cell_source_binning_v1.csv

import sys, os, io, atexit, gzip, warnings
import h5py
import numpy as np
import pandas as pd
import openpyxl
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

GENE_TABLE = os.environ.get(
    'PANEL68_XLSX',
    os.path.join(DATA_ROOT, 'data', 'processed',
                 'Supplementary_Table_NCC_68_genes.xlsx')
)
OUT_CSV = os.path.join(OUT_DIR, 'panel68_cell_source_binning_v1.csv')

COMPARTMENT_MAP = {
    'Tumor': 'Tumor',
    'Macro': 'Myeloid', 'Macro_MKI67': 'Myeloid', 'Mast': 'Myeloid',
    'Tcell': 'Lymphoid', 'Tcell_CD8': 'Lymphoid', 'Bcell': 'Lymphoid',
    'Plasma': 'Lymphoid',
    'Endo_PLVAP': 'Stromal', 'Endo_ACKR1': 'Stromal', 'Peri': 'Stromal',
    'vSMC': 'Stromal',
}
COMPARTMENTS = ['Tumor', 'Myeloid', 'Lymphoid', 'Stromal']
CELL_TYPES = list(COMPARTMENT_MAP.keys())


def merge_label(x):
    label_map = {
        '1:Tumor': 'Tumor', '2:Endo_PLVAP': 'Endo_PLVAP',
        '3:Endo_ACKR1': 'Endo_ACKR1', '4:Peri': 'Peri', '5:vSMC': 'vSMC',
        '6:Macro': 'Macro', '7:Macro_MKI67': 'Macro_MKI67', '8:Bcell': 'Bcell',
        '9:Tcell': 'Tcell', '10:Tcell_CD8': 'Tcell_CD8',
        '11:Tcell_CD8': 'Tcell_CD8', '12:Mast': 'Mast', '13:Plasma': 'Plasma',
    }
    return label_map.get(x, None)


def load_panel_genes():
    wb = openpyxl.load_workbook(GENE_TABLE)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    genes = []
    for r in rows[3:]:
        if r[0] is None:
            continue
        genes.append({'gene': str(r[0]).strip(),
                      'functional_layer': r[6],
                      'role': r[5]})
    df = pd.DataFrame(genes)
    print(f'Panel genes loaded: {len(df)}')
    return df


def main():
    panel = load_panel_genes()
    panel_genes = panel['gene'].tolist()

    anno = pd.read_csv(gzip.open(
        os.path.join(SCRNA_DIR, 'GSE159115_ccRCC_anno.csv.gz'), 'rt'))
    anno['celltype'] = anno['label'].apply(merge_label)
    anno = anno[anno['celltype'].notna()].copy()
    print(f'Annotated cells: {len(anno)}')

    accum = {g: {ct: [] for ct in CELL_TYPES} for g in panel_genes}
    found_genes = set()

    h5_files = sorted([f for f in os.listdir(SCRNA_DIR)
                       if f.startswith('GSM') and f.endswith('.h5')])
    for fn in h5_files:
        sample_id = fn.split('_SI_')[1].split('_')[0]
        sample_key = f'SI_{sample_id}'
        anno_sample = anno[anno['sample'] == sample_key]
        if len(anno_sample) == 0:
            print(f'  SKIP {sample_key}: no annotated cells')
            continue
        cell_to_ct = dict(zip(anno_sample['cell'], anno_sample['celltype']))
        cell_to_tot = dict(zip(anno_sample['cell'], anno_sample['total_UMI']))
        print(f'  Processing {sample_key}: {len(anno_sample)} annotated cells ...',
              flush=True)

        with h5py.File(os.path.join(SCRNA_DIR, fn), 'r') as f:
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

        gene_pos = {}
        for i, gn in enumerate(gene_names):
            if gn in accum and gn not in gene_pos:
                gene_pos[gn] = i
        found_genes.update(gene_pos.keys())

        cell_of_k = np.repeat(np.arange(n_cells), np.diff(indptr))
        mask = np.isin(indices, list(gene_pos.values()))
        cells_f = cell_of_k[mask]
        genes_f = indices[mask]
        data_f = data[mask]

        coo = sparse.coo_matrix(
            (data_f, (cells_f, genes_f)), shape=(n_cells, n_genes)).tocsr()

        ct_arr = np.array([cell_to_ct.get(f'{sample_key}_{bc}') for bc in barcodes],
                          dtype=object)
        tot_umi_arr = np.array(
            [cell_to_tot.get(f'{sample_key}_{bc}', 0.0) for bc in barcodes],
            dtype=float)

        for gn, gi in gene_pos.items():
            col = np.asarray(coo[:, gi].todense()).ravel().astype(np.float64)
            # per-cell CP10K normalization removes sequencing-depth confound
            norm = np.where(
                tot_umi_arr > 0,
                col / np.where(tot_umi_arr > 0, tot_umi_arr, 1.0) * 1e4,
                0.0)
            for ct in CELL_TYPES:
                ct_mask = ct_arr == ct
                accum[gn][ct].extend(norm[ct_mask].tolist())

    missing = [g for g in panel_genes if g not in found_genes]
    print(f'\nGenes found in h5 matrices: {len(found_genes)}/{len(panel_genes)}')
    if missing:
        print(f'  Missing: {missing}')

    ct_to_comp = COMPARTMENT_MAP
    rows = []
    for _, prow in panel.iterrows():
        gn = prow['gene']
        vals_by_comp = {}
        for comp in COMPARTMENTS:
            vals_by_comp[comp] = np.array(
                sum((accum[gn][ct] for ct in CELL_TYPES
                     if ct_to_comp[ct] == comp), []))
        tot = sum(float(v.sum()) for v in vals_by_comp.values())
        rec = {'gene': gn, 'functional_layer': prow['functional_layer'],
               'role': prow['role'], 'detected': tot > 0}
        means = {}
        for comp in COMPARTMENTS:
            vals = vals_by_comp[comp]
            n = len(vals)
            mean_umi = float(vals.mean()) if n > 0 else np.nan
            pct = float((vals > 0).mean() * 100) if n > 0 else np.nan
            share = float(vals.sum() / tot * 100) if tot > 0 else np.nan
            rec[f'{comp}_CP10K'] = round(mean_umi, 5) if not np.isnan(mean_umi) else ''
            rec[f'{comp}_pct_expr'] = round(pct, 2) if not np.isnan(pct) else ''
            rec[f'{comp}_CP10K_share_pct'] = round(share, 2) if not np.isnan(share) else ''
            means[comp] = mean_umi if n > 0 else 0.0
        valid_means = {k: v for k, v in means.items() if v > 0}
        if valid_means:
            dom = max(valid_means, key=valid_means.get)
            others = {k: v for k, v in valid_means.items() if k != dom}
            rec['dominant_compartment'] = dom
            rec['tumor_vs_max_other'] = round(
                means['Tumor'] / max(others.values()), 2) \
                if others and max(others.values()) > 0 else ''
        else:
            rec['dominant_compartment'] = 'not_detected'
            rec['tumor_vs_max_other'] = ''
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f'\nSaved: {OUT_CSV}')

    print('\nDominant compartment summary (detected genes):')
    det = out[out['detected']]
    print(det['dominant_compartment'].value_counts().to_string())
    print('\nPer-gene results:')
    show = ['gene', 'functional_layer', 'dominant_compartment',
            'Tumor_CP10K', 'Tumor_pct_expr', 'Myeloid_CP10K', 'Myeloid_pct_expr',
            'Lymphoid_CP10K', 'Lymphoid_pct_expr', 'Stromal_CP10K',
            'Stromal_pct_expr', 'tumor_vs_max_other']
    print(det[show].to_string(index=False))

    print('\nCompartment cell totals:')
    ct_counts = anno['celltype'].value_counts()
    comp_counts = {}
    for ct, n in ct_counts.items():
        comp_counts[COMPARTMENT_MAP[ct]] = comp_counts.get(COMPARTMENT_MAP[ct], 0) + n
    print(pd.Series(comp_counts).to_string())


if __name__ == '__main__':
    main()