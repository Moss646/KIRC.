#!/usr/bin/env python3
"""Run all Figure 2 clustering scripts in dependency order.

Each script is executed as a subprocess so its stdout appears in real time.
Any non-zero exit code stops the pipeline.

Usage:
    python run_upstream.py                  run every step
    python run_upstream.py --start 3        start from step 3
    python run_upstream.py --end 5          stop after step 5
    python run_upstream.py --only nmf_validation
    python run_upstream.py --list           list steps and exit
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

STEPS = [
    ("compute_clustering",         "consensus clustering (writes cache)"),
    ("build_classifier",           "train NCC classifier"),
    ("plot_figure2",               "Figure 2 main panels"),
    ("plot_figure2_pubsize",       "Figure 2 publication size"),
    ("supp_cdf_plot",              "supplementary CDF analysis"),
    ("supp_cdf_plot_pubsize",      "supplementary CDF publication size"),
    ("nmf_validation",             "NMF validation figure"),
    ("nmf_validation_metrics_v2",  "NMF vs KMeans metrics CSV"),
    ("nmf_validation_pubsize",     "NMF validation publication size"),
    ("table1_baseline",            "Table 1 baseline characteristics"),
    ("table_s2_ncc_classifier",    "Table S2 NCC classifier"),
    ("survival_blind_analysis_v1",
     "survival-blind consensus clustering + stability analysis"),
    ("survival_blind_supplementary_figure_v4",
     "survival-blind supplementary figure"),
]


def run_one(name, desc):
    script = HERE / f"{name}.py"
    if not script.exists():
        print(f"  [skip] {name}.py not found", flush=True)
        return 0

    print(f"\n>>> {name}.py  --  {desc}", flush=True)
    t0 = time.time()
    rc = subprocess.run([PY, str(script)], cwd=str(HERE)).returncode
    dt = time.time() - t0

    if rc != 0:
        print(f"    FAILED (exit {rc}, {dt:.0f}s)", flush=True)
    else:
        print(f"    OK ({dt:.0f}s)", flush=True)
    return rc


def main():
    parser = argparse.ArgumentParser(description="Figure 2 pipeline")
    parser.add_argument("--start", type=int, default=1,
                        help="first step (1-based)")
    parser.add_argument("--end", type=int, default=len(STEPS),
                        help="last step (inclusive)")
    parser.add_argument("--only", type=str, default=None,
                        help="run only one script by module name")
    parser.add_argument("--list", action="store_true",
                        help="list steps and exit")
    args = parser.parse_args()

    if args.list:
        for i, (name, desc) in enumerate(STEPS, start=1):
            print(f"{i:2d}  {name:<38s}  {desc}")
        return 0

    if args.only:
        match = [(n, d) for n, d in STEPS if n == args.only]
        if not match:
            print(f"unknown script: {args.only}")
            return 2
        name, desc = match[0]
        print(f"Figure 2 pipeline: running only {name}")
        return run_one(name, desc)

    start = max(1, args.start)
    end = min(len(STEPS), args.end)
    if start > end:
        print(f"invalid range: start={start}, end={end}")
        return 2

    print(f"Figure 2 pipeline: steps {start}-{end} of {len(STEPS)}")

    for i, (name, desc) in enumerate(STEPS, start=1):
        if i < start or i > end:
            continue
        print(f"\n{'=' * 64}")
        print(f"STEP {i}/{len(STEPS)}: {desc}")
        print(f"{'=' * 64}")

        rc = run_one(name, desc)
        if rc != 0:
            print(f"\nStopped at step {i} ({name}.py).")
            return rc

    print(f"\nAll steps {start}-{end} finished.")
    return 0


if __name__ == "__main__":
    sys.exit(main())