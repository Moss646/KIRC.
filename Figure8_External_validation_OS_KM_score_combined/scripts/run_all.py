"""Run the full NCC external validation pipeline.

Usage:
    python run.py                 run everything
    python run.py --skip-docx     skip typesetting
    python run.py --docx-only     typesetting only
    python run.py --only NAME     run only one script by filename
    python run.py --skip-prep     skip the mapping-JSON preparation step
"""

import os
import subprocess
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

from config import print_paths

PREP = [
    ("make_mapping_jsons.py",
     "prepare gene -> Ensembl mapping JSONs (needs Figure1 output)"),
]

CORE = [
    ("simple_cox_all_cohorts.py", "Simple Cox + log-rank, 3 cohorts"),
    ("prognostic_cox.py", "Stratified Cox + C-index"),
    ("external_validation_KM_pubsize.py", "KM curves"),
    ("classifier_concordance.py", "68-gene Pearson concordance"),
    ("key14_biological.py", "14-gene direction concordance + BH + Pearson"),
    ("panel_G_CPTAC_protein.py", "CPTAC protein validation"),
    ("diagnostics.py", "PH assumption + survival-time filter audit"),
]

EXTRA = [
    ("key14_cross_cohort_v2.py", "14-gene cross-cohort v2"),
    ("lipid16_cross_cohort_v2.py", "16-gene lipid remodeling cross-cohort v2"),
    ("lipid31_cross_cohort_v1.py", "30-gene combined cross-cohort v1"),
    ("nonncc9_fao_cross_cohort_v1.py", "9-gene non-NCC FAO + limma"),
    ("external_validation_KM_score_combined_v3.py", "Combined KM + FAO score"),
    ("driver_mut_fao_v1.py", "Driver mutations vs FAO score"),
]

SUMMARY = [
    ("verify_all_values.py", "Verify values against summary table"),
    ("generate_summary_csv.py", "Assemble summary CSV"),
]

DOCX = [
    (os.path.join("results", "tables", "build_table_docx.py"), "Main table docx"),
    (os.path.join("results", "tables", "build_split_tables.py"), "Split tables docx"),
    (os.path.join("results", "tables", "build_Table3.py"), "Reproducibility table docx"),
]

SUMMARY_CSV = os.path.join(
    REPO_ROOT, "results", "tables", "NCC_multi_cohort_validation_summary.csv"
)
CSV_GATE = "generate_summary_csv.py"


def run_step(python, label, abs_path, desc):
    print("\n" + "-" * 70)
    print(f"{label} {os.path.relpath(abs_path, REPO_ROOT)}")
    print(f"      {desc}")
    print("-" * 70)

    if not os.path.exists(abs_path):
        print(f"  SKIP: file not found: {abs_path}")
        return False

    t0 = time.time()
    try:
        result = subprocess.run(
            [python, abs_path],
            cwd=os.path.dirname(abs_path),
            capture_output=False,
            text=True,
        )
    except Exception as e:
        print(f"\n  EXCEPTION: {type(e).__name__}: {e}")
        return False

    elapsed = time.time() - t0

    if result.returncode != 0:
        print("\n" + "!" * 70)
        print(f"  FAILED: {os.path.basename(abs_path)}")
        print(f"  Exit code: {result.returncode}")
        print(f"  Elapsed  : {elapsed:.1f}s")
        print("!" * 70)
        return False

    print(f"\n  OK ({elapsed:.1f}s)")
    return True


def run_group(python, group, step_start, total_steps, failed):
    step = step_start
    for script, desc in group:
        step += 1
        abs_path = (
            os.path.join(REPO_ROOT, script)
            if script.startswith("results")
            else os.path.join(SCRIPT_DIR, script)
        )
        if not run_step(python, f"[{step}/{total_steps}]", abs_path, desc):
            failed.append(script)
    return step


def main():
    try:
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    except Exception:
        pass

    argv = sys.argv[1:]
    skip_docx = "--skip-docx" in argv
    docx_only = "--docx-only" in argv
    skip_prep = "--skip-prep" in argv

    only = None
    if "--only" in argv:
        i = argv.index("--only")
        if i + 1 < len(argv):
            only = argv[i + 1]

    print("=" * 70)
    print("NCC external validation - full pipeline")
    print("=" * 70)
    print()
    print_paths()
    print()

    if only:
        target = None
        for group in (PREP, CORE, EXTRA, SUMMARY):
            for script, desc in group:
                if script == only:
                    target = (script, desc)
                    break
            if target:
                break
        if not target:
            print(f"Unknown script: {only}")
            return 1
        ok = run_step(sys.executable, "[1/1]",
                      os.path.join(SCRIPT_DIR, target[0]), target[1])
        return 0 if ok else 1

    python = sys.executable
    failed = []
    start_all = time.time()

    prep = [] if (docx_only or skip_prep) else PREP
    core = [] if docx_only else CORE
    extra = [] if docx_only else EXTRA
    summary = [] if docx_only else SUMMARY
    docx = [] if skip_docx else DOCX

    total = len(prep) + len(core) + len(extra) + len(summary) + len(docx)
    print(f"Total steps: {total}")
    print()

    step = 0
    step = run_group(python, prep, step, total, failed)
    step = run_group(python, core, step, total, failed)
    step = run_group(python, extra, step, total, failed)
    step = run_group(python, summary, step, total, failed)

    if docx:
        if CSV_GATE in failed:
            print(f"\nSkipping docx: {CSV_GATE} failed.")
            docx = []
        elif not os.path.exists(SUMMARY_CSV):
            print(f"\nSkipping docx: summary CSV not found at {SUMMARY_CSV}")
            docx = []
            failed.append("(summary CSV missing)")

    step = run_group(python, docx, step, total, failed)

    elapsed = time.time() - start_all
    print("\n" + "=" * 70)
    print(f"Pipeline finished in {elapsed / 60:.1f} min")
    print(f"Steps attempted: {step}/{total}")
    if failed:
        print(f"Failed ({len(failed)}):")
        for f in failed:
            print(f"  - {f}")
        print("=" * 70)
        return 1
    print("All steps succeeded.")
    print(f"Results: {os.path.join(REPO_ROOT, 'results')}")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted by user.")
        sys.exit(130)
    except Exception:
        import traceback
        print("\n" + "!" * 70)
        print("Unhandled exception in run.py")
        print("!" * 70)
        traceback.print_exc()
        sys.exit(2)