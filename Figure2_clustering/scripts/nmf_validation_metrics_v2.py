# nmf_validation_metrics_v2.py
# Jaccard + consensus-matrix correlation metrics for NMF vs KMeans.

import sys, io, os, pickle, atexit
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import NMF
from sklearn.metrics import (adjusted_rand_score,
                             normalized_mutual_info_score,
                             cohen_kappa_score, jaccard_score)


def _wait_on_exit():
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.getch()
    except Exception:
        pass


atexit.register(_wait_on_exit)

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'data')
EXPR = os.path.join(DATA, 'KIRC_expr_log2_tpm.csv')
CACHE = os.path.join(REPO, 'cache', 'consensus_cache.pkl')
OUT = os.path.join(REPO, 'results',
                   'Supplementary_NMF_validation_metrics_v2.csv')

expr = pd.read_csv(EXPR, index_col=0)
sub = pd.read_csv(os.path.join(DATA, 'subtype_assignment_balanced.csv'))
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
with open(os.path.join(DATA, 'balanced_genes.txt')) as f:
    genes = [l.strip() for l in f if l.strip()]
genes = [g for g in genes if g in expr.index]
common = [s for s in expr.columns if s in sub_map]
X = expr.loc[genes, common].T.values
kmeans_labels = np.array([0 if sub_map[s] == 'Subtype_1' else 1
                          for s in common])
n_samples = len(common)
print(f'Data: {n_samples} samples x {len(genes)} genes')

best_err, best_labels, best_seed = 1e10, None, None
seed_labels = {}
for seed in range(100):
    nmf = NMF(n_components=2, init='nndsvda', max_iter=1000,
              random_state=seed)
    W = nmf.fit_transform(X)
    labels = np.argmax(W, axis=1)
    if len(set(labels)) < 2:
        continue
    seed_labels[seed] = labels.copy()
    if nmf.reconstruction_err_ < best_err:
        best_err = nmf.reconstruction_err_
        best_labels = labels.copy()
        best_seed = seed
print(f'Best seed: {best_seed}, err={best_err:.4f}')

n0 = (best_labels == 0).sum()
n1 = (best_labels == 1).sum()
flip = (n0 > n1)
nmf_labels = (1 - best_labels) if flip else best_labels


def align_to_best(lab):
    cm = np.zeros((2, 2), dtype=int)
    for a, b in zip(lab, best_labels):
        cm[a, b] += 1
    r, c = linear_sum_assignment(-cm)
    m = {int(a): int(b) for a, b in zip(r, c)}
    return np.array([m[int(v)] for v in lab])


co = np.zeros((n_samples, n_samples))
for lab in seed_labels.values():
    al = align_to_best(lab)
    co += (al[:, None] == al[None, :]).astype(float)
nmf_consensus = co / len(seed_labels)

ari = adjusted_rand_score(kmeans_labels, nmf_labels)
nmi = normalized_mutual_info_score(kmeans_labels, nmf_labels)
kappa = cohen_kappa_score(kmeans_labels, nmf_labels)
agr = (kmeans_labels == nmf_labels).mean()
j_per = jaccard_score(kmeans_labels, nmf_labels, average=None)
j_macro = jaccard_score(kmeans_labels, nmf_labels, average='macro')

cm = np.zeros((2, 2), dtype=int)
for kl, nl in zip(kmeans_labels, nmf_labels):
    cm[kl, nl] += 1

s1 = nmf_labels == 0
s2 = nmf_labels == 1
w1 = nmf_consensus[np.ix_(s1, s1)][np.triu_indices(s1.sum(), 1)].mean()
w2 = nmf_consensus[np.ix_(s2, s2)][np.triu_indices(s2.sum(), 1)].mean()
bt = nmf_consensus[np.ix_(s1, s2)].mean()

with open(CACHE, 'rb') as f:
    d = pickle.load(f)
kk = 2 if 2 in d['results'] else d['best_k']
km_cons = d['results'][kk]['consensus']
cache_samples = d['samples']
sidx = {s: i for i, s in enumerate(cache_samples)}
cidx = [sidx[s] for s in common]
km_cons = km_cons[np.ix_(cidx, cidx)]
print(f'Cache best_k={d["best_k"]}; using k={kk} consensus matrix '
      f'({km_cons.shape[0]} samples aligned)')

iu = np.triu_indices(n_samples, k=1)
a = nmf_consensus[iu]
b = km_cons[iu]
r_pe = pearsonr(a, b)
r_sp = spearmanr(a, b)

rows = [
    ('n_samples', n_samples),
    ('n_genes', len(genes)),
    ('best_seed', best_seed),
    ('reconstruction_error', round(best_err, 4)),
    ('NMF_S1_size', int((nmf_labels == 0).sum())),
    ('NMF_S2_size', int((nmf_labels == 1).sum())),
    ('ARI', round(ari, 4)),
    ('NMI', round(nmi, 4)),
    ('Cohens_kappa', round(kappa, 4)),
    ('Overall_agreement', round(agr, 4)),
    ('Agreement_count', f'{int(agr * n_samples)}/{n_samples}'),
    ('Confusion_S1S1', cm[0, 0]),
    ('Confusion_S1toS2', cm[0, 1]),
    ('Confusion_S2toS1', cm[1, 0]),
    ('Confusion_S2S2', cm[1, 1]),
    ('Jaccard_S1', round(float(j_per[0]), 4)),
    ('Jaccard_S2', round(float(j_per[1]), 4)),
    ('Jaccard_macro', round(float(j_macro), 4)),
    ('NMF_consensus_within_S1_mean', round(float(w1), 4)),
    ('NMF_consensus_within_S2_mean', round(float(w2), 4)),
    ('NMF_consensus_between_mean', round(float(bt), 4)),
    ('Consensus_corr_Pearson_r', round(float(r_pe[0]), 4)),
    ('Consensus_corr_Pearson_p', float(f'{r_pe[1]:.3e}')),
    ('Consensus_corr_Spearman_rho', round(float(r_sp[0]), 4)),
    ('Consensus_corr_Spearman_p', float(f'{r_sp[1]:.3e}')),
]

os.makedirs(os.path.dirname(OUT), exist_ok=True)
pd.DataFrame(rows, columns=['metric', 'value']).to_csv(OUT, index=False)

print('\n=== RESULTS ===')
for k, v in rows:
    print(f'{k}: {v}')
print(f'\nSaved: {OUT}')