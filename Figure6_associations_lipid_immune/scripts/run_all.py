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

# config.py resolves the expression matrix as
# MEL_DATA_ROOT/tcga_kirc/KIRC_expr_log2_tpm.csv, so point MEL_DATA_ROOT
# at this project's data/ directory.
os.environ['MEL_DATA_ROOT'] = os.path.join(root, 'data')
os.environ['NOPAUSE'] = '1'

script_order = [
    'metab_immune_corr_v4.py',
    'fig_Figure6_main_v1.py',
]

search_dirs = [base, os.path.join(root, 'scripts')]

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

for name in script_order:
    path = find_script(name)
    print(f"Running {name}...")
    result = subprocess.run([python, path], cwd=root)
    if result.returncode != 0:
        print(f"Failed at: {name}")
        sys.exit(result.returncode)

print("All done. Figures in results/")