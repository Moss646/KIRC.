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

script_names = [
    'prepare_oncopredict_inputs.py',
    'run_oncoPredict_calcPheno.R',
    'run_tidepy.py',
    'Fig7_main_pubsize_v2.py',
    'Fig7_PanelA_violin_faceted.py',
]

search_dirs = [
    base,
    os.path.join(base, 'scripts'),
    os.path.join(base, 'scripts', 'R'),
    os.path.join(root, 'scripts'),
]

def find_script(name):
    for d in search_dirs:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None

missing = [s for s in script_names if find_script(s) is None]
if missing:
    print("Missing script(s):")
    for m in missing:
        print(" ", m)
    sys.exit(1)

python = sys.executable
rscript = 'Rscript'

steps = [
    ("prepare inputs", [python, find_script('prepare_oncopredict_inputs.py')]),
    ("oncoPredict (R)", [rscript, find_script('run_oncoPredict_calcPheno.R')]),
    ("TIDE analysis", [python, find_script('run_tidepy.py')]),
    ("Figure 7 main", [python, find_script('Fig7_main_pubsize_v2.py')]),
    ("Panel A violin", [python, find_script('Fig7_PanelA_violin_faceted.py')]),
]

for name, cmd in steps:
    print(f"Running {name}...")
    result = subprocess.run(cmd, cwd=root)
    if result.returncode != 0:
        print(f"Failed at: {name}")
        sys.exit(result.returncode)

print("All done. Figures in results/figures/")