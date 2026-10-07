#!/usr/bin/env python3
import os
import sys
import subprocess

base = os.path.dirname(os.path.abspath(__file__))
root = base
while not os.path.isdir(os.path.join(root, 'data')):
    parent = os.path.dirname(root)
    if parent == root:
        break
    root = parent

os.environ['MEL_DATA_ROOT'] = root

script_order = [
    'extract_spill.R',
    '03_xCell_analysis.R',
    'run_cibersort_official.R',
    'run_ssgsea_immune.R',
    'make_purity_strata_v3.py',
    'make_panel68_xlsx.py',
    '01_pseudobulk_ncc_v5.py',
    '02_celltype_subtype_expression_v5.py',
    'panel68_cell_source_binning_v1.py',
    'hif_gene_sets_zscore_v1.py',
    'run_ssgsea_hif_sets_v1.R',
    'run_gsva_hif_sets_v1.R',
    'run_camera_hif_sets_v1.R',
    'hif_gene_sets_per_gene_check_v1.py',
    'hif_gene_sets_figure_v1.py',
    'hif_pathway_activity_v2.py',
    'run_immune_genes_strata_limma_v3.R',
    'make_stratified_mwu_v2.py',
    'fig4_immune_merged_v17.py',
    'immune_features_composition_independence_v4.py',
    'lipid_remodeling_boxplot_v5.py',
    'pla2_inflammatory_boxplot_v1.py',
    'run_xcell_v2_pubsize.py',
    '03_plot_v5_pubsize.py',
]

search_dirs = [
    base,
    os.path.join(base, 'scripts'),
    os.path.join(root, 'scripts'),
]

def find_script(name):
    for d in search_dirs:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None

missing = [s for s in script_order if find_script(s) is None]
if missing:
    print("Missing script(s):")
    for m in missing:
        print(" ", m)
    sys.exit(1)

python = sys.executable
rscript = 'Rscript'

for name in script_order:
    path = find_script(name)
    if name.endswith('.R'):
        cmd = [rscript, path]
    else:
        cmd = [python, path]
    print(f"Running {name}...")
    result = subprocess.run(cmd, cwd=root)
    if result.returncode != 0:
        print(f"Failed at: {name}")
        sys.exit(result.returncode)

print("All done. Figures in results/figures/")