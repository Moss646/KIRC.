#!/usr/bin/env python3
"""Univariate Cox screening for lipid metabolism genes."""

import io
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DATA_DIR = REPO_ROOT / "data"

EXPR_FILE = DATA_DIR / "tcga_kirc" / "KIRC_expr_log2_tpm.csv"
CLIN_FILE = DATA_DIR / "tcga_kirc" / "KIRC_clinical_cleaned.csv"
GENE_POOL_FILE = DATA_DIR / "msigdb_lipid_gene_pool_dedup.csv"
OUTPUT_DIR = DATA_DIR / "cox_screening"
P_VALUE_THRESHOLD = 0.05


def cox_univariate(x, time, status):
    mask = ~(np.isnan(x) | np.isnan(time) | np.isnan(status))
    x = x[mask]
    time = time[mask]
    status = status[mask].astype(int)

    if len(x) < 10 or int(status.sum()) < 3:
        return np.nan, np.nan, np.nan, np.nan, False

    sx = x.std()
    if sx < 1e-10:
        return np.nan, np.nan, np.nan, np.nan, False

    xz = (x - x.mean()) / sx
    df = pd.DataFrame({"time": time, "status": status, "expr": xz})

    try:
        cph = CoxPHFitter()
        cph.fit(df, duration_col="time", event_col="status")
        beta = float(cph.params_["expr"])
        se = float(cph.summary.loc["expr", "se(coef)"])
        p_value = float(cph.summary.loc["expr", "p"])
        hr = float(np.exp(beta))
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, False

    return beta, se, hr, p_value, True


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty(n)
    out[order] = q
    return out


def main():
    print("Univariate Cox screening")

    if not EXPR_FILE.exists():
        raise SystemExit(f"Expression file not found: {EXPR_FILE}")
    if not CLIN_FILE.exists():
        raise SystemExit(f"Clinical file not found: {CLIN_FILE}")

    expr_df = pd.read_csv(EXPR_FILE, index_col=0)
    clinical = pd.read_csv(CLIN_FILE, index_col=0)
    gene_pool = pd.read_csv(GENE_POOL_FILE)
    lipid_genes_all = gene_pool["Gene_Symbol"].tolist()

    common_genes = [g for g in lipid_genes_all if g in expr_df.index]
    print(f"expression matrix: {expr_df.shape[0]} genes x {expr_df.shape[1]} samples")
    print(f"lipid genes in expression matrix: {len(common_genes)}/{len(lipid_genes_all)}")

    clin_valid = clinical.dropna(subset=["os", "os_time_days"])
    clin_valid = clin_valid[clin_valid["os_time_days"] > 0]

    common_samples = sorted(set(expr_df.columns) & set(clin_valid.index))
    os_time = clin_valid.loc[common_samples, "os_time_days"].values.astype(float)
    os_status = clin_valid.loc[common_samples, "os"].values.astype(int)

    print(f"samples: {len(common_samples)}, events: {int(os_status.sum())}")

    expr_sub = expr_df.loc[common_genes, common_samples].values.astype(float)
    results = []
    t0 = time.time()

    for i, gene in enumerate(common_genes):
        if (i + 1) % 500 == 0:
            print(f"{i + 1}/{len(common_genes)} ({time.time() - t0:.1f}s)")

        beta, se, hr, pval, ok = cox_univariate(
            expr_sub[i, :].flatten(), os_time, os_status
        )
        if not ok:
            continue

        results.append({
            "gene": gene,
            "coef": beta,
            "HR": hr,
            "se": se,
            "p_value": pval,
            "significant": int(pval < P_VALUE_THRESHOLD),
        })

    if not results:
        raise RuntimeError("no gene converged")

    results_df = pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)
    results_df["FDR"] = bh_fdr(results_df["p_value"].values)
    results_df["significant_fdr"] = (results_df["FDR"] < 0.05).astype(int)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results_df.to_csv(OUTPUT_DIR / "cox_results_all.csv", index=False)
    results_df[results_df["significant"] == 1].to_csv(
        OUTPUT_DIR / "cox_results_significant.csv", index=False
    )
    results_df[results_df["significant_fdr"] == 1].to_csv(
        OUTPUT_DIR / "cox_results_fdr_significant.csv", index=False
    )

    with open(OUTPUT_DIR / "cox_sig_genes_list.txt", "w") as f:
        for g in results_df.loc[results_df["significant"] == 1, "gene"]:
            f.write(g + "\n")

    summary = {
        "total_tested": len(common_genes),
        "total_converged": len(results_df),
        "p_significant": int(results_df["significant"].sum()),
        "fdr_significant": int(results_df["significant_fdr"].sum()),
        "n_risk": int((results_df["HR"] > 1).sum()),
        "n_protective": int((results_df["HR"] < 1).sum()),
    }
    with open(OUTPUT_DIR / "cox_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"saved {len(results_df)} genes to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()