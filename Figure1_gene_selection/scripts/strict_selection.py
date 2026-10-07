#!/usr/bin/env python3
"""Strict multi-criteria selection of lipid metabolism genes."""

import io
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DATA_DIR = REPO_ROOT / "data"

COX_ALL = DATA_DIR / "cox_screening" / "cox_results_all.csv"
EXPR_FILE = DATA_DIR / "tcga_kirc" / "KIRC_expr_log2_tpm.csv"
GENE_POOL_DEDUP = DATA_DIR / "msigdb_lipid_gene_pool_dedup.csv"
GENESETS_DEDUP = DATA_DIR / "msigdb_lipid_genesets_dedup.csv"
OUTPUT_DIR = DATA_DIR / "gene_selection"

P_THRESHOLD = 0.001
LOG2HR_THRESHOLD = 0.30

LAYER_KEYWORDS = {
    "Core_Lipid_Enzyme": [
        "FATTY_ACID_METABOLISM", "FATTY_ACID_BIOSYNTHETIC", "FATTY_ACID_BETA",
        "LIPID_BIOSYNTHETIC", "LIPID_CATABOLIC", "CHOLESTEROL_BIOSYNTHESIS",
        "PHOSPHOLIPID", "TRIGLYCERIDE", "GLYCEROLIPID", "SPHINGOLIPID",
        "STEROID_BIOSYNTHESIS", "BILE_ACID", "ACYL_COA", "FATTY_ACID_ELONGATION",
        "CARNITINE", "KETONE_BODY", "LIPASE", "PEROXISOME",
    ],
    "Lipid_Droplet_Storage": [
        "LIPID_DROPLET", "LIPID_STORAGE", "LIPID_LOCALIZATION",
        "LIPID_TRANSPORT", "LIPID_HOMEOSTASIS",
    ],
    "Transcriptional_Regulator": [
        "ADIPOGENESIS", "PPAR", "SREBF", "SREBP", "LXR",
        "STEROID_HORMONE_RECEPTOR", "REGULATION_OF_LIPID",
    ],
    "Ferroptosis_ROS": [
        "FERROPTOSIS", "LIPID_OXIDATION", "LIPID_PEROXIDATION",
        "REACTIVE_OXYGEN", "GLUTATHIONE",
    ],
    "Signaling_Lipoprotein": [
        "LIPOPROTEIN", "CHYLOMICRON", "HDL", "LDL", "VLDL",
        "CHOLESTEROL_EFFLUX", "PHOSPHOLIPID_TRANSPORT",
        "PROTEIN_LIPID_COMPLEX", "PHOSPHATIDYLINOSITOL",
        "PROSTAGLANDIN", "EICOSANOID", "ARACHIDONIC",
        "PHOSPHOLIPASE", "LIPID_SIGNALING",
    ],
}


def main():
    if not EXPR_FILE.exists():
        raise SystemExit(f"Expression file not found: {EXPR_FILE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cox = pd.read_csv(COX_ALL)
    for col in ["p_value", "FDR", "se"]:
        if col in cox.columns:
            cox[col] = pd.to_numeric(
                cox[col].astype(str).str.replace(r"[\[\]]", "", regex=True),
                errors="coerce",
            )
    for col in ["HR", "coef"]:
        if col in cox.columns:
            cox[col] = pd.to_numeric(cox[col], errors="coerce")

    cox = cox.dropna(subset=["p_value", "HR"])

    expr = pd.read_csv(EXPR_FILE, index_col=0)
    gene_pool = pd.read_csv(GENE_POOL_DEDUP)
    lipid_genes_in_expr = [g for g in gene_pool["Gene_Symbol"].tolist() if g in expr.index]

    print(f"cox results: {len(cox)}")
    print(f"lipid genes in expression matrix: {len(lipid_genes_in_expr)}")

    df = cox[cox["gene"].isin(lipid_genes_in_expr)].copy()
    print(f"start: {len(df)}")

    df["log2HR"] = np.log2(df["HR"].values)

    before = len(df)
    df = df[df["p_value"] < P_THRESHOLD]
    print(f"p < {P_THRESHOLD}: {before} -> {len(df)}")

    before = len(df)
    df = df[df["log2HR"].abs() > LOG2HR_THRESHOLD]
    print(f"|log2HR| > {LOG2HR_THRESHOLD}: {before} -> {len(df)}")

    before = len(df)
    surviving_genes = df["gene"].tolist()
    expr_sub = expr.loc[surviving_genes].values
    gene_sd = np.std(expr_sub, axis=1, ddof=1)
    median_sd = np.median(gene_sd)
    df["sd"] = gene_sd
    df = df[df["sd"] > median_sd]
    print(f"SD > median ({median_sd:.3f}): {before} -> {len(df)}")

    before = len(df)
    selected_genes = df["gene"].tolist()
    expr_sel = expr.loc[selected_genes].values.T
    corr_matrix = np.corrcoef(expr_sel, rowvar=False)
    n_genes = len(selected_genes)
    removed = set()

    for i in range(n_genes):
        if i in removed:
            continue
        for j in range(i + 1, n_genes):
            if j in removed:
                continue
            if abs(corr_matrix[i, j]) > 0.85:
                pi = df.iloc[i]["p_value"]
                pj = df.iloc[j]["p_value"]
                if pi <= pj:
                    removed.add(j)
                else:
                    removed.add(i)
                    break

    keep_mask = [i for i in range(len(df)) if i not in removed]
    df = df.iloc[keep_mask].reset_index(drop=True)
    print(f"|r| < 0.85: {before} -> {len(df)}")

    genesets = pd.read_csv(GENESETS_DEDUP)
    selected_set = set(df["gene"].tolist())
    gene_layer = {}

    for gene in selected_set:
        gene_sets = genesets[genesets["Genes"].str.contains(rf"\b{gene}\b", na=False)]
        layer_scores = {layer: 0 for layer in LAYER_KEYWORDS}

        for _, gs in gene_sets.iterrows():
            gs_name = gs["Geneset_Name"].upper()
            for layer, keywords in LAYER_KEYWORDS.items():
                if any(kw.upper() in gs_name for kw in keywords):
                    layer_scores[layer] += 1

        best_layer = max(layer_scores, key=layer_scores.get)
        gene_layer[gene] = best_layer if layer_scores[best_layer] > 0 else "Other"

    layer_counts = Counter(gene_layer[g] for g in selected_set)

    out_cols = ["gene", "HR", "log2HR", "p_value", "FDR", "coef", "se", "sd"]
    out_df = df[out_cols].copy()
    out_df["functional_layer"] = out_df["gene"].map(gene_layer)
    out_df.to_csv(OUTPUT_DIR / "selected_genes_strict.csv", index=False)

    with open(OUTPUT_DIR / "selected_genes_list.txt", "w") as f:
        for g in df["gene"].tolist():
            f.write(g + "\n")

    nmf_expr = expr.loc[df["gene"].tolist()]
    nmf_expr.to_csv(OUTPUT_DIR / "nmf_input_expr.csv")

    log_data = {
        "input_genes": len(lipid_genes_in_expr),
        "final_selected": len(df),
        "p_threshold": P_THRESHOLD,
        "log2hr_threshold": LOG2HR_THRESHOLD,
        "n_risk": int((df["HR"] > 1).sum()),
        "n_protective": int((df["HR"] < 1).sum()),
        "layer_distribution": dict(layer_counts),
    }
    with open(OUTPUT_DIR / "selection_log.json", "w") as f:
        json.dump(log_data, f, indent=2)

    print(f"selected {len(df)} genes")
    print(f"output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()