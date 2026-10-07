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
    '02_differential_exp.R',
    '01_convert_ensembl_to_symbol.R',
    'make_msigdb_gmt.py',
    '03_fgsea.R',
    'gsva_R.R',
    'gsva_limma.R',
    'Fig3_main_v18.py',
    'cptac_lipid_remodeling_check_v2.py',
]

search_dirs = [
    base,
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
    cmd = [rscript, path] if name.endswith('.R') else [python, path]
    print(f"Running {name}...")
    result = subprocess.run(cmd, cwd=root)
    if result.returncode != 0:
        print(f"Failed at: {name}")
        sys.exit(result.returncode)

print("All done. Figures in results/figures/")