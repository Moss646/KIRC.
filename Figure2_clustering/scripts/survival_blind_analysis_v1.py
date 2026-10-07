# -*- coding: utf-8 -*-
"""
survival_blind_analysis_v1.py
=============================
Merged upstream analysis for the survival-blind clustering (reviewer
response): combines survival_blind_clustering_v1.py (main k = 2 experiment)
and survival_blind_k_stability_v3_pacfix.py (k = 2-6 stability) into a
single script. Both original scripts and their outputs are left untouched;
the merged script reproduces their numbers exactly (same seed 20261002 and
identical subsampling logic, so the k = 2 partition matches v1 sample by
sample).

Design (survival-blind by construction)
---------------------------------------
* Gene universe: full MSigDB lipid gene pool after IC dedup, intersected
  with the expressed TCGA-KIRC matrix. No survival information enters the
  clustering at any step.
* Simplified Monti consensus clustering: for each k in 2-6, 1,000 iterations;
  each iteration subsamples 80% of samples and 80% of genes, re-z-scores
  within the subsample, and runs k-means. Pairwise co-clustering frequencies
  form the consensus matrix; final partition = average linkage on 1 - C.
* Per-k stability: PAC (all pairs; and within-cluster pairs for reference),
  mean silhouette on the consensus distance (1 - C) and in expression space
  (Euclidean, z-scored genes), ARI/NMI vs the published S1/S2 assignment,
  chi-square test.
* k = 2 main experiment: majority-vote label alignment, Cohen's kappa,
  concordance rate, direct-k-means sensitivity check, and OS log-rank of the
  blind clusters.

Outputs
-------
- results/survival_blind_analysis_v1_stats.csv          (per-k stability table)
- results/survival_blind_analysis_v1_k2_concordance.csv (k = 2 metrics)
- results/survival_blind_analysis_v1_assignment.csv     (per-sample k = 2 labels)
- results/survival_blind_analysis_v1_crosstab.csv       (2 x 2 table)
- results/survival_blind_analysis_v1_consensus_k{k}.npy (cached matrices)
"""

import io
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans
from sklearn.metrics import (adjusted_rand_score, cohen_kappa_score,
                             normalized_mutual_info_score, silhouette_score)
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(HERE)

EXPR_FILE = os.path.join(PROJ_ROOT, "data", "KIRC_expr_log2_tpm.csv")
CLIN_FILE = os.path.join(PROJ_ROOT, "data", "KIRC_clinical_cleaned.csv")
POOL_FILE = os.path.join(PROJ_ROOT, "data",
                         "msigdb_lipid_gene_pool_dedup.csv")
SUB_FILE = os.path.join(PROJ_ROOT, "data",
                        "subtype_assignment_balanced.csv")
OUT_DIR = os.path.join(PROJ_ROOT, "results")

K_RANGE = [2, 3, 4, 5, 6]
N_ITER = 1000
SAMPLE_FRAC = 0.80
GENE_FRAC = 0.80
SEED = 20261002
PAC_LO, PAC_HI = 0.1, 0.9


def consensus_matrix(X, k, n_iter, sample_frac, gene_frac, seed):
    """Simplified Monti consensus clustering. X: samples x genes (z-scored)."""
    n_samp, n_gene = X.shape
    ns = int(round(sample_frac * n_samp))
    ng = int(round(gene_frac * n_gene))
    co = np.zeros((n_samp, n_samp), dtype=np.int64)
    cnt = np.zeros((n_samp, n_samp), dtype=np.int64)
    rng = np.random.RandomState(seed)
    for it in range(n_iter):
        s_idx = rng.choice(n_samp, ns, replace=False)
        g_idx = rng.choice(n_gene, ng, replace=False)
        sub = X[np.ix_(s_idx, g_idx)]
        sub = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-10)
        km = KMeans(n_clusters=k, n_init=3, random_state=rng.randint(1 << 31))
        lab = km.fit_predict(sub)
        idx = np.ix_(s_idx, s_idx)
        co[idx] += (lab[:, None] == lab[None, :]).astype(np.int64)
        cnt[idx] += 1
        if (it + 1) % 250 == 0:
            print(f"    k={k}: iteration {it + 1}/{n_iter}", flush=True)
    return co / np.maximum(cnt, 1)


def pac_values(C, lab):
    """PAC over all upper-triangle pairs and over within-cluster pairs only."""
    iu = np.triu_indices_from(C, k=1)
    vals = C[iu]
    amb = (vals > PAC_LO) & (vals < PAC_HI)
    pac_all = float(amb.mean())
    same = lab[:, None] == lab[None, :]
    within = same[iu]
    denom = int(within.sum())
    pac_within = float(amb[within].mean()) if denom > 0 else np.nan
    return pac_all, pac_within, denom


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- inputs ------------------------------------------------------------
    pool = pd.read_csv(POOL_FILE)["Gene_Symbol"].astype(str).tolist()
    expr = pd.read_csv(EXPR_FILE, index_col=0)
    genes_ok = [g for g in pool if g in expr.index]
    print(f"gene pool {len(pool)} -> expressed: {len(genes_ok)}", flush=True)
    Xz = expr.loc[genes_ok].T
    Xz = (Xz - Xz.mean(axis=0)) / (Xz.std(axis=0) + 1e-10)
    sub = pd.read_csv(SUB_FILE).set_index("Patient")["Subtype"]
    clin = pd.read_csv(CLIN_FILE).set_index("patient_id")
    common = [s for s in Xz.index if s in sub.index and s in clin.index]
    Xz = Xz.loc[common]
    samples = np.array(common)
    y_sub = sub.loc[samples].to_numpy()
    time = clin.loc[samples, "os_time_days"].to_numpy(float)
    event = clin.loc[samples, "os"].to_numpy(float).astype(int)
    n_samp = len(samples)
    print(f"samples: {n_samp} (S1={(y_sub == 'Subtype_1').sum()}, "
          f"S2={(y_sub == 'Subtype_2').sum()})", flush=True)

    Xmat = Xz.to_numpy(float)

    # ---- per-k consensus clustering + stability metrics --------------------
    rows = []
    k2_labels_raw = None
    for k in K_RANGE:
        npy = os.path.join(
            OUT_DIR, f"survival_blind_analysis_v1_consensus_k{k}.npy")
        if os.path.exists(npy):
            print(f"consensus matrix k = {k}: loading cached npy ...",
                  flush=True)
            C = np.load(npy)
        else:
            print(f"consensus clustering k = {k} ...", flush=True)
            C = consensus_matrix(Xmat, k, N_ITER, SAMPLE_FRAC, GENE_FRAC, SEED)
            np.save(npy, C)

        D_cond = squareform(np.clip(1.0 - C, 0, 1), checks=False)
        lab = fcluster(linkage(D_cond, method="average"), t=k,
                       criterion="maxclust")
        pac_all, pac_within, within_pairs = pac_values(C, lab)
        sil_cons = float(silhouette_score(np.clip(1.0 - C, 0, 1), lab,
                                          metric="precomputed"))
        sil_expr = float(silhouette_score(Xmat, lab, metric="euclidean"))
        ari = adjusted_rand_score(y_sub, lab)
        nmi = normalized_mutual_info_score(y_sub, lab)
        ct_k = pd.crosstab(pd.Series(y_sub, name="Published"),
                           pd.Series(lab, name="Blind"))
        chi2, chi_p, _, _ = chi2_contingency(ct_k.to_numpy())
        rows.append({
            "k": k,
            "PAC_allpairs": round(pac_all, 4),
            "PAC_within": round(pac_within, 4),
            "within_pairs": within_pairs,
            "mean_silhouette_consensus_distance": round(sil_cons, 4),
            "mean_silhouette_expression": round(sil_expr, 4),
            "ARI_vs_published": round(float(ari), 4),
            "NMI_vs_published": round(float(nmi), 4),
            "chi2": round(float(chi2), 3),
            "chi2_p": chi_p,
        })
        if k == 2:
            k2_labels_raw = lab
    stats = pd.DataFrame(rows)
    stats.to_csv(os.path.join(
        OUT_DIR, "survival_blind_analysis_v1_stats.csv"), index=False)
    print(stats.to_string(index=False), flush=True)

    # ---- k = 2 main experiment ---------------------------------------------
    lab_raw = k2_labels_raw
    align = {}
    for c in np.unique(lab_raw):
        mask = lab_raw == c
        vals, counts = np.unique(y_sub[mask], return_counts=True)
        align[c] = vals[np.argmax(counts)]
    clus = np.array(["S1-blind" if align[c] == "Subtype_1" else "S2-blind"
                     for c in lab_raw])
    m1 = clus == "S1-blind"
    m2 = clus == "S2-blind"

    ari2 = adjusted_rand_score(y_sub, clus)
    clus_mapped = np.where(m1, "Subtype_1", "Subtype_2")
    kappa = cohen_kappa_score(y_sub, clus_mapped)
    ct = pd.crosstab(pd.Series(y_sub, name="Published"),
                     pd.Series(clus, name="Survival-blind"))
    chi2_2, chi_p_2, _, _ = chi2_contingency(ct.to_numpy())
    n_conc = int(ct.to_numpy()[0, 0] + ct.to_numpy()[1, 1])
    conc = n_conc / n_samp

    km_full = KMeans(n_clusters=2, n_init=20,
                     random_state=SEED).fit_predict(Xmat)
    ari_km = adjusted_rand_score(y_sub, km_full)
    km_al = np.array([None] * n_samp, dtype=object)
    for c in np.unique(km_full):
        mask = km_full == c
        vals, counts = np.unique(y_sub[mask], return_counts=True)
        km_al[mask] = vals[np.argmax(counts)]
    kappa_km = cohen_kappa_score(y_sub, km_al)

    lr = logrank_test(time[m1], time[m2], event[m1], event[m2])

    k2 = pd.DataFrame([
        ("gene_pool_size", len(pool)),
        ("genes_expressed", len(genes_ok)),
        ("n_samples", n_samp),
        ("n_iter", N_ITER),
        ("sample_fraction", SAMPLE_FRAC),
        ("gene_fraction", GENE_FRAC),
        ("ARI_vs_published", round(float(ari2), 4)),
        ("kappa_vs_published", round(float(kappa), 4)),
        ("concordance", round(conc, 4)),
        ("concordant_n", n_conc),
        ("chi2_crosstab", round(float(chi2_2), 3)),
        ("chi2_p", chi_p_2),
        ("ARI_direct_kmeans", round(float(ari_km), 4)),
        ("kappa_direct_kmeans", round(float(kappa_km), 4)),
        ("blind_S1_n", int(m1.sum())),
        ("blind_S2_n", int(m2.sum())),
        ("blind_S1_median_OS_days", float(np.median(time[m1]))),
        ("blind_S2_median_OS_days", float(np.median(time[m2]))),
        ("logrank_p_OS", lr.p_value),
    ], columns=["metric", "value"])
    k2.to_csv(os.path.join(
        OUT_DIR, "survival_blind_analysis_v1_k2_concordance.csv"), index=False)
    ct.to_csv(os.path.join(
        OUT_DIR, "survival_blind_analysis_v1_crosstab.csv"))
    pd.DataFrame({
        "Patient": samples, "Published": y_sub,
        "Survival_blind_cluster": clus,
        "OS_time_days": time, "OS_event": event,
    }).to_csv(os.path.join(
        OUT_DIR, "survival_blind_analysis_v1_assignment.csv"), index=False)
    print(ct.to_string())
    print(k2.to_string(index=False), flush=True)
    print("merged analysis complete", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())