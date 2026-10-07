# nmf_validation.py
# NMF validation of the KMeans consensus clustering subtypes.

import sys, os, io, atexit
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.decomposition import NMF
from sklearn.metrics import (adjusted_rand_score,
                             normalized_mutual_info_score,
                             cohen_kappa_score)


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

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['svg.fonttype'] = 'none'

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'data')
FIG = os.path.join(REPO, 'results')
os.makedirs(FIG, exist_ok=True)

EXPR = os.path.join(DATA, 'KIRC_expr_log2_tpm.csv')

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

print(f'Data: {len(common)} samples x {len(genes)} genes')
print(f'KMeans: S1={(kmeans_labels == 0).sum()}, '
      f'S2={(kmeans_labels == 1).sum()}')

best_err, best_H = 1e10, None
best_labels, best_seed = None, None
for seed in range(100):
    nmf = NMF(n_components=2, init='nndsvda', max_iter=1000,
              random_state=seed)
    W = nmf.fit_transform(X)
    labels = np.argmax(W, axis=1)
    if len(set(labels)) < 2:
        continue
    if nmf.reconstruction_err_ < best_err:
        best_err = nmf.reconstruction_err_
        best_H = nmf.components_.copy()
        best_labels = labels.copy()
        best_seed = seed

print(f'Best seed: {best_seed}, reconstruction error: {best_err:.4f}')

# smaller cluster = S1
n0 = (best_labels == 0).sum()
n1 = (best_labels == 1).sum()
flip = (n0 > n1)
nmf_labels = (1 - best_labels) if flip else best_labels
s1_comp, s2_comp = (1, 0) if flip else (0, 1)
H = np.array([best_H[s1_comp], best_H[s2_comp]])

ari = adjusted_rand_score(kmeans_labels, nmf_labels)
nmi = normalized_mutual_info_score(kmeans_labels, nmf_labels)
kappa = cohen_kappa_score(kmeans_labels, nmf_labels)
agr = (kmeans_labels == nmf_labels).sum() / len(nmf_labels)
cm = np.zeros((2, 2), dtype=int)
for kl, nl in zip(kmeans_labels, nmf_labels):
    cm[kl, nl] += 1

print(f'\n=== RESULTS ===')
print(f'NMF: S1={np.sum(nmf_labels == 0)}, '
      f'S2={np.sum(nmf_labels == 1)}')
print(f'Agreement: {agr:.1%} ({int(agr * len(common))}/{len(common)})')
print(f'ARI: {ari:.4f}, NMI: {nmi:.4f}')
print(f"Cohen's kappa: {kappa:.4f}")
print(f'S1 stability: {cm[0,0]}/{cm[0,0]+cm[0,1]} '
      f'= {cm[0,0]/(cm[0,0]+cm[0,1])*100:.1f}%')
print(f'S2 stability: {cm[1,1]}/{cm[1,0]+cm[1,1]} '
      f'= {cm[1,1]/(cm[1,0]+cm[1,1])*100:.1f}%')

C1, C2 = '#C65A5A', '#4C78A8'
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
plt.subplots_adjust(wspace=0.35)
fig.patch.set_facecolor('white')

# A: NMF component difference
ax = axes[0]
diff = H[0] - H[1]
order = np.argsort(diff)
ax.barh(range(len(genes)), diff[order],
        color=[C1 if v > 0 else C2 for v in diff[order]],
        alpha=0.85, height=0.7)
ax.axvline(0, color='black', lw=0.8)
ax.set_yticks([])
ax.set_ylabel(f'NCC genes (n = {len(genes)})',
              fontsize=28, fontweight='bold')
ax.set_xlabel('NMF component difference (S1 - S2)',
              fontsize=28, fontweight='bold')
ax.set_title('A', fontsize=32, fontweight='bold', loc='left', pad=6)
for spine in ax.spines.values():
    spine.set_linewidth(2)
ax.tick_params(labelsize=18)
ax.legend(handles=[Patch(facecolor=C1, label='S1 (NMF)'),
                   Patch(facecolor=C2, label='S2 (NMF)')],
          loc='lower right', fontsize=20, framealpha=0.95,
          edgecolor='black')

# B: confusion matrix
ax = axes[1]
im = ax.imshow(cm, cmap='Blues')
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(['S1 (NMF)', 'S2 (NMF)'],
                   fontsize=28, fontweight='bold')
ax.set_yticklabels(['S1 (KMeans)', 'S2 (KMeans)'],
                   fontsize=28, fontweight='bold')
ax.set_xlabel('NMF-K-means subtype concordance',
              fontsize=28, fontweight='bold')
for i in range(2):
    for j in range(2):
        ax.text(j, i, str(cm[i, j]), ha='center', va='center',
                fontsize=28, fontweight='bold',
                color='white' if cm[i, j] > cm.max() / 2 else 'black')
ax.set_title('B', fontsize=32, fontweight='bold', loc='left', pad=6)
for spine in ax.spines.values():
    spine.set_linewidth(2)
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Number of samples', fontsize=20, fontweight='bold')

plt.tight_layout()
fig.savefig(os.path.join(FIG, 'Supplementary_NMF_validation_v3.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, 'Supplementary_NMF_validation_v3.svg'),
            format='svg', bbox_inches='tight', facecolor='white')
plt.close()
print(f'\nSaved: Supplementary_NMF_validation_v3.png/.svg')