#!/usr/bin/env python3
"""Run the Figure 1 upstream pipeline."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PY = sys.executable

STEPS = [
    ("download_msigdb_lipid", "MSigDB lipid gene-set download"),
    ("go_dedup_lipid", "GO:BP hierarchical dedup"),
    ("cox_screening", "univariate Cox screening"),
    ("strict_selection", "multi-criteria strict selection"),
    ("balance_genes_for_clustering", "68-gene balanced panel"),
]


def run(script_name):
    cmd = [PY, str(HERE / f"{script_name}.py")]
    return subprocess.run(cmd).returncode


def verify_classifier(classifier_path):
    print("STEP 6: regenerate figure and verify classifier")
    rc = run("gene_selection_fig_v12_pubsize")
    if rc != 0:
        print(f"figure generation failed with code {rc}")
        sys.exit(rc)

    strict = pd.read_csv(REPO / "data/gene_selection/selected_genes_strict.csv")
    balanced = pd.read_csv(
        REPO / "data/consensus_cluster/balanced_genes.txt",
        header=None,
    ).iloc[:, 0].tolist()

    if not classifier_path:
        print("OFFICIAL_CLASSIFIER not set; skipping classifier verification.")
        print("strict selected genes:", len(strict))
        print("balanced genes:", len(balanced))
        return

    classifier_path = Path(classifier_path)
    if not classifier_path.exists():
        print(f"Classifier file not found: {classifier_path}; skipping verification.")
        return

    with open(classifier_path) as f:
        off = set(json.load(f)["genes"])

    print("strict selected genes:", len(strict))
    print("balanced genes:", len(balanced))
    print("official classifier:", len(off))
    print("official subset of balanced:", off.issubset(set(balanced)))
    print("official subset of strict:", off.issubset(set(strict["gene"])))

    if off == set(balanced):
        print("68-gene classifier matches official set")
    else:
        print("68-gene classifier differs from official set")
        print("only official:", sorted(off - set(balanced))[:20])
        print("only new:", sorted(set(balanced) - off)[:20])


def main():
    parser = argparse.ArgumentParser(description="Figure 1 upstream pipeline")
    parser.add_argument("--start", type=int, default=1, help="start step 1..5")
    parser.add_argument("--end", type=int, default=5, help="end step 1..5")
    parser.add_argument(
        "--no-figure",
        dest="figure",
        action="store_false",
        help="skip figure generation and classifier verification",
    )
    parser.add_argument(
        "--classifier",
        default=os.environ.get("OFFICIAL_CLASSIFIER"),
        help="Path to official subtype_classifier.json (optional)",
    )
    args = parser.parse_args()

    for i, (mod_name, desc) in enumerate(STEPS, start=1):
        if i < args.start or i > args.end:
            continue
        print(f"STEP {i}/5: {desc}")
        rc = run(mod_name)
        if rc != 0:
            print(f"step {i} ({mod_name}) failed with code {rc}")
            sys.exit(rc)

    print(f"pipeline finished: steps {args.start}-{args.end}")

    if args.figure and args.end >= 5:
        verify_classifier(args.classifier)
    elif args.figure and args.end < 5:
        print("figure skipped because --end < 5")


if __name__ == "__main__":
    main()