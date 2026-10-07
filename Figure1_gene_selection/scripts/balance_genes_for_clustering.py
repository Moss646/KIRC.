#!/usr/bin/env python3
"""Select a balanced 68-gene panel from strict-selected genes."""

import os
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DATA_DIR = REPO_ROOT / "data"

SEL_FILE = DATA_DIR / "gene_selection" / "selected_genes_strict.csv"
OUT_DIR = DATA_DIR / "consensus_cluster"
OUT_FILE = OUT_DIR / "balanced_genes.txt"
OLD_FILE = OUT_DIR / "balanced_genes_prev.txt"

QUOTA = {
    "Core_Lipid_Enzyme": 20,
    "Signaling_Lipoprotein": 15,
    "Transcriptional_Regulator": 15,
    "Lipid_Droplet_Storage": 5,
    "Ferroptosis_ROS": 2,
    "Other": 11,
}
TARGET = sum(QUOTA.values())
assert TARGET == 68, f"quota sum is {TARGET}, expected 68"


def main():
    sel = pd.read_csv(SEL_FILE)
    selected = []

    for layer, n in QUOTA.items():
        sub = sel[sel["functional_layer"] == layer].sort_values("p_value")
        take = min(n, len(sub))
        selected.extend(sub.head(take)["gene"].tolist())
        print(f"{layer}: {len(sub)} available, took {take}/{n}")

    if len(selected) < TARGET:
        picked = set(selected)
        remaining = sel[~sel["gene"].isin(picked)].sort_values("p_value")
        need = TARGET - len(selected)
        extra = remaining.head(need)["gene"].tolist()
        selected.extend(extra)
        print(f"added {need} genes: {extra}")

    selected = selected[:TARGET]

    if OLD_FILE.exists():
        with open(OLD_FILE) as f:
            old = {line.strip() for line in f if line.strip()}
        new = set(selected)
        print(f"overlap with previous: {len(old & new)}/{TARGET}")
        if old - new:
            print("only previous:", sorted(old - new))
        if new - old:
            print("only new:", sorted(new - old))

    gene_layer = dict(zip(sel["gene"], sel["functional_layer"]))
    dist = Counter(gene_layer[g] for g in selected)
    print("layer distribution:")
    for layer in QUOTA:
        print(f"  {layer}: {dist.get(layer, 0)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_FILE, "w") as f:
        for g in selected:
            f.write(g + "\n")

    print(f"wrote {OUT_FILE} ({len(selected)} genes)")


if __name__ == "__main__":
    main()