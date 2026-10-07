"""Assemble NCC_multi_cohort_validation_summary.csv from provenance JSON sidecars.

Final step of run_all.py. Each analysis script writes a JSON sidecar into
results/tables/_provenance/<script>.json via config.emit_provenance. This
collector merges them into the summary table and adds per-cell provenance
columns (source script, data version, timestamp).

The previous hand-maintained table, if present, is used for validation only.
"""

import csv
import json
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", "tables"))
PROVENANCE_DIR = os.path.join(TABLES_DIR, "_provenance")
MANUAL_CSV = os.path.join(TABLES_DIR, "NCC_multi_cohort_validation_summary.manual.csv")
OUT_CSV = os.path.join(TABLES_DIR, "NCC_multi_cohort_validation_summary.csv")

COHORTS = ["TCGA-KIRC (Training)", "E-MTAB-1980", "ICGC RECA-EU", "CPTAC"]

ROWS = [
    ("Cohort characteristics", "S1 / S2"),
    ("Cohort characteristics", "S1 proportion"),
    ("Prognostic validation", "NCC HR"),
    ("Prognostic validation", "Cox P"),
    ("Prognostic validation", "Log-rank P"),
    ("Prognostic validation", "Stratified NCC HR (Stage+Grade+Age)"),
    ("Prognostic validation", "Stratified NCC 95% CI"),
    ("Prognostic validation", "Stratified NCC P"),
    ("Prognostic validation", "N in stratified model (complete Stage+Grade+Age)"),
    ("Prognostic validation", "Stratified C-index (Harrell's C)"),
    ("Prognostic validation", "C-index 95% CI (bootstrap 2000)"),
    ("Classifier robustness", "NCC genes available"),
    ("Classifier robustness", "Missing genes"),
    ("Classifier robustness", "68-gene Pearson r vs TCGA"),
    ("Classifier robustness", "68-gene 95% CI"),
    ("Classifier robustness", "68-gene Pearson P"),
    ("Biological validation", "14-gene direction concordance"),
    ("Biological validation", "14-gene BH-sig genes (MWU FDR<0.05)"),
    ("Biological validation", "14-gene Pearson r vs TCGA"),
    ("Biological validation", "14-gene 95% CI"),
    ("Biological validation", "14-gene Pearson P"),
]


def _norm(v):
    if v is None:
        return "-"
    v = str(v).strip()
    # strip leading zeros from scientific-notation exponents: 1.8e-04 -> 1.8e-4
    v = re.sub(r"e([+-])0+(\d)", r"e\1\2", v)
    # normalize ranges to en-dash: 1.39-2.67 -> 1.39–2.67
    # (minus signs in exponents like 1.8e-4 are left intact)
    v = re.sub(r"(\d)-(\d)", lambda m: m.group(1) + "\u2013" + m.group(2), v)
    return v


def _load():
    cells = {}
    if not os.path.isdir(PROVENANCE_DIR):
        print(f"[COLLECTOR] no provenance dir: {PROVENANCE_DIR}")
        return cells

    for fn in sorted(os.listdir(PROVENANCE_DIR)):
        if not fn.endswith(".json"):
            continue
        try:
            with open(os.path.join(PROVENANCE_DIR, fn), encoding="utf-8") as f:
                d = json.load(f)
        except Exception as e:
            print(f"[COLLECTOR] skip {fn}: {e}")
            continue

        for r in d.get("rows", []):
            key = (r["section"], r["item"], r["cohort"])
            cells[key] = {
                "value": _norm(r.get("value")),
                "script": d.get("script"),
                "data_version": d.get("data_version"),
                "generated_at": d.get("generated_at"),
            }

    return cells


def _main():
    cells = _load()

    s12 = {}
    for coh in COHORTS:
        c = cells.get(("Cohort characteristics", "S1 / S2", coh))
        if not c:
            continue
        m = re.match(r"(\d+)\s*of\s*(\d+)", c["value"])
        if m:
            s1, s2 = int(m.group(1)), int(m.group(2))
            s12[coh] = (s1, s2, c["script"])

    out_rows = []
    for section, item in ROWS:
        row = {"Section": section, "Item": item}
        src = {}
        dv = set()
        gen = set()

        for coh in COHORTS:
            if item == "S1 proportion":
                if coh in s12:
                    s1, s2, sc = s12[coh]
                    row[coh] = f"{100.0 * s1 / (s1 + s2):.1f}%"
                    src[coh] = f"derived from S1/S2 ({sc})"
                else:
                    row[coh] = "-"
                    src[coh] = "-"
            else:
                c = cells.get((section, item, coh))
                if c:
                    row[coh] = c["value"]
                    src[coh] = c["script"]
                    if c["data_version"]:
                        dv.add(c["data_version"])
                    if c["generated_at"]:
                        gen.add(c["generated_at"])
                else:
                    row[coh] = "-"
                    src[coh] = "-"

        row["Src:TCGA"] = src.get("TCGA-KIRC (Training)", "-")
        row["Src:E-MTAB"] = src.get("E-MTAB-1980", "-")
        row["Src:ICGC"] = src.get("ICGC RECA-EU", "-")
        row["Src:CPTAC"] = src.get("CPTAC", "-")
        row["Data version"] = sorted(dv)[0] if dv else "unknown"
        row["Generated"] = sorted(gen)[0] if gen else "unknown"
        out_rows.append(row)

    cols = [
        "Section",
        "Item",
        *COHORTS,
        "Src:TCGA",
        "Src:E-MTAB",
        "Src:ICGC",
        "Src:CPTAC",
        "Data version",
        "Generated",
    ]

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print(f"[COLLECTOR] wrote {len(out_rows)} rows -> {OUT_CSV}")

    if os.path.isfile(MANUAL_CSV):
        _validate(OUT_CSV, MANUAL_CSV)


def _cell_equal(a, b):
    a, b = _norm(a), _norm(b)
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) < 1e-6
    except ValueError:
        return False


def _validate(auto, manual):
    def read_csv(p):
        with open(p, encoding="utf-8") as f:
            return list(csv.reader(f))

    a = read_csv(auto)
    m = read_csv(manual)
    mism = []

    for i in range(1, min(len(a), len(m))):
        for j in range(6):
            av = a[i][j] if j < len(a[i]) else ""
            mv = m[i][j] if j < len(m[i]) else ""
            if not _cell_equal(av, mv):
                mism.append((i, j, mv, av))

    if mism:
        print(f"[COLLECTOR] VALIDATION: {len(mism)} data cell(s) differ from manual:")
        for i, j, mv, av in mism:
            print(f"  row{i + 1} col{j}: manual='{mv}'  auto='{av}'")
    else:
        print("[COLLECTOR] VALIDATION: all 21 data rows match manual CSV")


if __name__ == "__main__":
    _main()