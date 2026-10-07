# compute_clustering.py
# Recompute the consensus-clustering cache from TCGA-KIRC.
# On the first run, when subtype_assignment_balanced.csv does not yet exist,
# it is written from the clustering result (smaller cluster = Subtype_1).

import atexit
import io
import os
import pickle
import re
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, cophenet, fcluster
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")


def _wait_on_exit():
    try:
        if os.name == "nt":
            import msvcrt
            msvcrt.getch()
    except Exception:
        pass


atexit.register(_wait_on_exit)

if sys.platform == "win32" and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
CACHE = os.path.join(REPO, "cache", "consensus_cache.pkl")
DATA_DIR = os.path.join(REPO, "data")
np.random.seed(42)


def _pick(df, names):
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    for n in names:
        for lc, c in lower.items():
            if n.lower() in lc:
                return c
    return None


def _to_patient(x):
    s = str(x).strip()
    m = re.match(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", s)
    return m.group(1) if m else s[:12]


def _norm_stage(x):
    s = str(x).strip().upper()
    m = re.search(r"\b(IV|III|II|I)\b", s)
    return m.group(1) if m else None


def _norm_os(x):
    s = str(x).strip().lower()
    if s in ("1", "true", "dead", "deceased"):
        return 1
    if s in ("0", "false", "alive", "living"):
        return 0
    try:
        return int(float(s))
    except Exception:
        return None


def _load_clinical():
    info_path = os.path.join(DATA_DIR, "TCGA-KIRC_clinical_info.csv")
    surv_path = os.path.join(DATA_DIR, "TCGA-KIRC_clinical_survival.csv")
    for p in (info_path, surv_path):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    info = pd.read_csv(info_path)
    surv = pd.read_csv(surv_path)

    iid = _pick(info, ["bcr_patient_barcode", "submitter_id",
                       "case_id", "patient"])
    istage = _pick(info, ["ajcc_pathologic_tumor_stage",
                          "ajcc_pathologic_stage", "tumor_stage", "stage"])
    sid = _pick(surv, ["bcr_patient_barcode", "submitter_id",
                       "case_id", "sample"])
    sos = _pick(surv, ["OS", "vital_status"])
    stime = _pick(surv, ["OS.time", "os_time_days",
                         "days_to_death", "days_to_last_follow_up"])

    missing = [n for n, v in [
        ("info.id", iid), ("info.stage", istage),
        ("surv.id", sid), ("surv.OS", sos), ("surv.OS.time", stime),
    ] if v is None]

    if missing:
        print("\n[column match failed] missing:", missing)
        print("info columns:", list(info.columns))
        print("surv columns:", list(surv.columns))
        raise SystemExit(1)

    info2 = pd.DataFrame({
        "Patient": info[iid].map(_to_patient),
        "stage":   info[istage].map(_norm_stage),
    }).drop_duplicates("Patient").set_index("Patient")

    surv2 = pd.DataFrame({
        "Patient":      surv[sid].map(_to_patient),
        "os":           surv[sos].map(_norm_os),
        "os_time_days": pd.to_numeric(surv[stime], errors="coerce"),
    }).drop_duplicates("Patient").set_index("Patient")

    merged = info2.join(surv2, how="inner").dropna(
        subset=["os", "os_time_days"])
    print(f"  Clinical loaded: {len(merged)} patients "
          f"(stage missing in {merged['stage'].isna().sum()})")
    return merged


def _align_to_reference(labels, samples):
    """Relabel clusters to match the reference subtype file
    (Subtype_1 -> 0, Subtype_2 -> 1)."""
    ref_path = os.path.join(DATA_DIR, "subtype_assignment_balanced.csv")
    ref = pd.read_csv(ref_path)
    ref_map = dict(zip(ref["Patient"], ref["Subtype"]))

    ref_labels = np.array([
        0 if ref_map.get(s, "Subtype_2") == "Subtype_1" else 1
        for s in samples
    ])

    n_clusters = int(labels.max()) + 1
    cm = np.zeros((n_clusters, 2), dtype=int)
    for a, b in zip(labels, ref_labels):
        cm[int(a), int(b)] += 1

    r, c = linear_sum_assignment(-cm)
    remap = {int(rr): int(cc) for rr, cc in zip(r, c)}
    return np.array([remap[int(v)] for v in labels])


def _write_reference_from_labels(labels, samples):
    """First run: label clusters, smaller cluster = Subtype_1, write file."""
    n0 = int((labels == 0).sum())
    n1 = int((labels == 1).sum())
    if n0 > n1:
        labels = 1 - labels

    ref_path = os.path.join(DATA_DIR, "subtype_assignment_balanced.csv")
    out = pd.DataFrame({
        "Patient": samples,
        "Subtype": ["Subtype_1" if v == 0 else "Subtype_2" for v in labels],
    })
    out.to_csv(ref_path, index=False)
    n_s1 = int((labels == 0).sum())
    n_s2 = int((labels == 1).sum())
    print(f"  First run: wrote {ref_path} (S1={n_s1}, S2={n_s2})")


def run_clustering():
    expr_path = os.path.join(DATA_DIR, "KIRC_expr_log2_tpm.csv")
    if not os.path.exists(expr_path):
        raise FileNotFoundError(expr_path)

    EXPR = pd.read_csv(expr_path, index_col=0)
    CLIN = _load_clinical()

    with open(os.path.join(REPO, "data", "balanced_genes.txt")) as f:
        balanced_genes = [line.strip() for line in f if line.strip()]
    sel_genes = [g for g in balanced_genes if g in EXPR.index]
    print(f"  Genes: {len(sel_genes)}")

    X = EXPR.loc[sel_genes].T.values
    X = StandardScaler().fit_transform(X)
    samples = list(EXPR.columns)
    n_samples = X.shape[0]

    results = {}
    for k in range(2, 7):
        t1 = time.time()
        co = np.zeros((n_samples, n_samples))
        cnt = np.zeros((n_samples, n_samples))
        for run in range(100):
            n_sub = max(int(n_samples * 0.8), k * 5)
            idx = np.random.choice(n_samples, n_sub, replace=False)
            try:
                km = KMeans(n_clusters=k, init="k-means++", n_init=3,
                            max_iter=300, random_state=run)
                ls = km.fit_predict(X[idx, :])
                for i in range(n_sub):
                    for j in range(i, n_sub):
                        cnt[idx[i], idx[j]] += 1
                        cnt[idx[j], idx[i]] += 1
                        if ls[i] == ls[j]:
                            co[idx[i], idx[j]] += 1
                            co[idx[j], idx[i]] += 1
            except Exception:
                pass

        consensus = np.divide(co, cnt, out=np.zeros_like(co), where=cnt > 0)
        dist = 1.0 - consensus
        np.fill_diagonal(dist, 0)
        Z = linkage(squareform(dist, checks=False), method="average")
        labels = fcluster(Z, k, criterion="maxclust") - 1

        try:
            sil = float(silhouette_score(X, labels)
                        if len(set(labels)) > 1 else 0)
        except Exception:
            sil = 0
        pac = float(np.mean((consensus > 0.1) & (consensus < 0.9)))
        coph = float(cophenet(Z, squareform(dist, checks=False))[0])
        score = sil + coph - pac

        print(f"    k={k}: sil={sil:.3f}, coph={coph:.3f}, "
              f"PAC={pac:.3f} -> {score:.3f} ({time.time()-t1:.0f}s)")
        results[k] = {"consensus": consensus, "labels": labels,
                      "silhouette": sil, "cophenetic": coph,
                      "pac": pac, "score": score}

    best_k = max(results,
                 key=lambda k: results[k]["score"]
                 if len(set(results[k]["labels"])) >= 2 else -1)
    labels = results[best_k]["labels"]
    print(f"  Best k={best_k}")

    ref_path = os.path.join(DATA_DIR, "subtype_assignment_balanced.csv")
    if os.path.exists(ref_path):
        labels = _align_to_reference(labels, samples)
    else:
        _write_reference_from_labels(labels, samples)

    results[best_k]["labels"] = labels

    idx_map = {s: i for i, s in enumerate(samples)}
    os_rows = []
    for pid, row in CLIN.iterrows():
        if pid in idx_map:
            i = idx_map[pid]
            os_rows.append({"patient": pid, "subtype": labels[i],
                            "os": row["os"],
                            "os_time": row["os_time_days"],
                            "stage": row["stage"]})
    os_df = pd.DataFrame(os_rows)
    for s in sorted(set(labels)):
        n = (labels == s).sum()
        print(f"    S{s+1}: {n} patients")

    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump({"results": results, "best_k": best_k, "labels": labels,
                     "os_df": os_df, "samples": samples}, f)
    print(f"  Cached: {CACHE}")
    return results, best_k, labels, os_df


print("=" * 50)
print("  Consensus clustering - force regeneration")
print("=" * 50)

if os.path.exists(CACHE):
    os.remove(CACHE)
    print(f"  Removed old cache: {CACHE}")

results, best_k, labels, os_df = run_clustering()
print("  Done. Cache written; no figure generated by this script.")