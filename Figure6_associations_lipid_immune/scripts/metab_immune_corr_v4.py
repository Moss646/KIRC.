# metab_immune_corr_v4.py
# Lipid-metabolism x immune-feature correlation analysis.
# 30 genes (14 FA metabolism + 16 phospholipid-remodeling) x 7 immune
# features = 210 pairs. Spearman, partial Spearman (ImmuneScore +
# StromaScore), subtype-stratified Spearman. BH across all 210 within each
# scheme. Output: results/metabolic_immune_corr_30genes_v1.csv

import sys, os, io, atexit
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, pearsonr
from scipy.stats import t as tdist

from config import EXPR, SUBTYPE_CSV, TIDE_CSV, XCELL_CSV, require_large_inputs
import config


def _hold():
    if os.environ.get('NOPAUSE'):
        return
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.getch()
        else:
            input('Press Enter to exit...')
    except (EOFError, KeyboardInterrupt):
        pass


atexit.register(_hold)

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

OUT = os.path.join(config.RESULTS_DIR, "metabolic_immune_corr_30genes_v1.csv")

PANEL_14 = [
    ("CPT1A", "14g fatty acid metabolism", "FAO", "S1 lower"),
    ("CPT2", "14g fatty acid metabolism", "FAO", "S1 lower"),
    ("ACOX1", "14g fatty acid metabolism", "FAO", "S1 lower"),
    ("ACADM", "14g fatty acid metabolism", "FAO", "S1 lower"),
    ("ACADSB", "14g fatty acid metabolism", "FAO", "S1 lower"),
    ("ACAA2", "14g fatty acid metabolism", "FAO", "S1 lower"),
    ("SLC27A2", "14g fatty acid metabolism", "uptake", "S1 lower"),
    ("CD36", "14g fatty acid metabolism", "uptake", "S1 lower"),
    ("PLIN2", "14g fatty acid metabolism", "storage", "S1 lower"),
    ("PPARG", "14g fatty acid metabolism", "regulator", "S1 lower"),
    ("PPARGC1A", "14g fatty acid metabolism", "regulator", "S1 lower"),
    ("FASN", "14g fatty acid metabolism", "synthesis", "S1 higher (modest)"),
    ("SCD", "14g fatty acid metabolism", "synthesis", "ns"),
    ("ITPKA", "14g fatty acid metabolism", "signaling", "S1 higher"),
]
PANEL_16 = [
    ("LPCAT1", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("LPCAT2", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("LPCAT3", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("LPCAT4", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("MBOAT1", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("MBOAT2", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("LCLAT1", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("MBOAT7", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("PLA2G4A", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("PLA2G6", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("PNPLA8", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("PLA2G4C", "16g phospholipid remodeling", "tumor-intrinsic enzyme"),
    ("PLA2G2A", "16g phospholipid remodeling", "inflammatory PLA2"),
    ("PLA2G1B", "16g phospholipid remodeling", "inflammatory PLA2"),
    ("PLA2G2D", "16g phospholipid remodeling", "inflammatory PLA2"),
    ("PLA2G4F", "16g phospholipid remodeling", "inflammatory PLA2"),
]
PANEL = ([(g, p, grp, "") for g, p, grp in PANEL_16]
         + [(g, p, grp, d) for g, p, grp, d in PANEL_14])
GENES = [g for g, _, _, _ in PANEL]
CKPT = ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD274"]
TIDE_COL = "TIDE"


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    ranked_fdr = np.minimum.accumulate(
        ranked[::-1] * n / np.arange(n, 0, -1))[::-1]
    ranked_fdr = np.minimum(ranked_fdr, 1.0)
    out = np.empty(n)
    out[order] = ranked_fdr
    return out


def partial_spearman(x, y, covs):
    """Rank-based partial correlation: residualize ranks on covariate ranks."""
    rx = rankdata(np.asarray(x, dtype=float))
    ry = rankdata(np.asarray(y, dtype=float))
    rc = [rankdata(np.asarray(c, dtype=float)) for c in covs]
    A = np.column_stack([np.ones_like(rx)] + rc).astype(np.float64)
    k = A.shape[1] - 1

    def resid(v):
        v = np.asarray(v, dtype=np.float64)
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    ex, ey = resid(rx), resid(ry)
    rho, _ = pearsonr(ex, ey)
    n = len(x)
    if abs(rho) >= 1.0:
        return rho, 0.0
    t = rho * np.sqrt((n - 2 - k) / (1 - rho ** 2))
    p = 2 * tdist.sf(abs(t), n - 2 - k)
    return rho, p


def main():
    require_large_inputs()

    expr = pd.read_csv(EXPR).set_index("gene").T
    expr.index.name = "sample"
    expr.columns = [c.upper() for c in expr.columns]
    expr = expr.reset_index()
    missing = [g for g in GENES if g not in expr.columns]
    assert not missing, f"genes missing from expression matrix: {missing}"

    sub = pd.read_csv(SUBTYPE_CSV).rename(columns={"Patient": "sample"})
    sub["subtype"] = sub["Subtype"].map({"Subtype_1": "S1",
                                         "Subtype_2": "S2"})
    tide = pd.read_csv(TIDE_CSV).rename(columns={"Unnamed: 0": "sample"})
    xcell = pd.read_csv(XCELL_CSV).rename(columns={"Unnamed: 0": "celltype"})

    imm = xcell[xcell["celltype"] == "ImmuneScore"].iloc[0, 1:].rename("ImmuneScore")
    imm = imm.reset_index().rename(columns={"index": "sample"})
    strm = xcell[xcell["celltype"] == "StromaScore"].iloc[0, 1:].rename("StromaScore")
    strm = strm.reset_index().rename(columns={"index": "sample"})

    df = (sub[["sample", "subtype"]]
          .merge(tide[["sample", TIDE_COL]], on="sample", how="inner")
          .merge(expr[["sample"] + GENES + CKPT], on="sample", how="inner")
          .merge(imm, on="sample", how="inner")
          .merge(strm, on="sample", how="inner"))
    assert len(df) == 533 and (df["subtype"] == "S1").sum() == 172
    print("merged n =", len(df), "| S1 =", (df["subtype"] == "S1").sum())

    # Coerce all numeric columns to float, drop any rows with non-numeric
    # entries that survived the merges (TIDE / xCell sometimes carry
    # empty-string or "NA" cells from their source exports).
    num_cols = (GENES + CKPT + ["ImmuneScore", "StromaScore", TIDE_COL])
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    dropped = df[num_cols].isna().any(axis=1).sum()
    if dropped:
        print(f"dropping {dropped} rows with non-numeric values in "
              f"expression / immune columns")
        df = df.dropna(subset=num_cols).reset_index(drop=True)
    print("analysis n =", len(df), "| S1 =", (df["subtype"] == "S1").sum())

    df["purity"] = 1.0 - (df["ImmuneScore"] + df["StromaScore"])

    rows = []
    for gene, panel, grp, s1dir in PANEL:
        for immvar in [TIDE_COL] + CKPT:
            rho_all, p_all = spearmanr(df[gene], df[immvar])
            d1 = df[df["subtype"] == "S1"]
            d2 = df[df["subtype"] == "S2"]
            rho_s1, p_s1 = spearmanr(d1[gene], d1[immvar])
            rho_s2, p_s2 = spearmanr(d2[gene], d2[immvar])
            rho_adj1, p_adj1 = partial_spearman(
                df[gene], df[immvar], [df["ImmuneScore"]])
            rho_adj2, p_adj2 = partial_spearman(
                df[gene], df[immvar], [df["ImmuneScore"], df["StromaScore"]])
            rows.append({
                "panel": panel, "gene_group": grp, "metabolic_gene": gene,
                "S1_direction": s1dir,
                "immune_variable": immvar,
                "n_all": len(df), "rho_all": rho_all, "p_all": p_all,
                "n_S1": len(d1), "rho_S1": rho_s1, "p_S1": p_s1,
                "n_S2": len(d2), "rho_S2": rho_s2, "p_S2": p_s2,
                "rho_adjImmune": rho_adj1, "p_adjImmune": p_adj1,
                "rho_adjPurity": rho_adj2, "p_adjPurity": p_adj2,
            })
    res = pd.DataFrame(rows)
    res["fdr_all"] = bh_fdr(res["p_all"].values)
    res["fdr_adjImmune"] = bh_fdr(res["p_adjImmune"].values)
    res["fdr_adjPurity"] = bh_fdr(res["p_adjPurity"].values)

    order = ["panel", "gene_group", "metabolic_gene", "S1_direction",
             "immune_variable", "n_all", "rho_all", "p_all", "fdr_all",
             "n_S1", "rho_S1", "p_S1", "rho_S2", "p_S2",
             "rho_adjImmune", "p_adjImmune", "fdr_adjImmune",
             "rho_adjPurity", "p_adjPurity", "fdr_adjPurity"]
    fmt = res[order].copy()
    for c in fmt.columns:
        if c not in ("panel", "gene_group", "metabolic_gene", "S1_direction",
                     "immune_variable", "n_all", "n_S1", "n_S2"):
            fmt[c] = fmt[c].map(lambda v: f"{v:.4g}")
    fmt.to_csv(OUT, index=False)
    print("saved:", OUT)

    npairs = len(res)
    sig_all = set(res[res["fdr_all"] < 0.05].index)
    sig_new = set(res[res["fdr_adjPurity"] < 0.05].index)
    print("\nFDR<0.05: full=%d/%d | adjPurity(Imm+Stroma)=%d/%d"
          % (len(sig_all), npairs, len(sig_new), npairs))
    for panel in sorted(res["panel"].unique()):
        m = res["panel"] == panel
        print("%s: full=%d/%d | adjPurity=%d/%d"
              % (panel,
                 (m & (res["fdr_all"] < 0.05)).sum(), m.sum(),
                 (m & (res["fdr_adjPurity"] < 0.05)).sum(), m.sum()))

    m14 = res["panel"] == "14g fatty acid metabolism"
    print("\n14-gene subset within 210-family: adjPurity sig = %d/%d"
          % ((m14 & (res["fdr_adjPurity"] < 0.05)).sum(), m14.sum()))
    tide14 = res[m14 & (res["immune_variable"] == "TIDE")]
    print("TIDE x 14 genes, adjPurity FDR<0.05: %d/14"
          % (tide14["fdr_adjPurity"] < 0.05).sum())


if __name__ == '__main__':
    main()