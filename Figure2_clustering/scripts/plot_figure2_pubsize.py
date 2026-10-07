# plot_figure2_pubsize.py
# Publication-size version of Figure 2 (180 x 135 mm @ 600 dpi).
# Same analysis as plot_figure2.py; only the figure block differs.

import sys, os, io, atexit, warnings, pickle
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')


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
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
FIG = os.path.join(REPO, 'results')
CACHE = os.path.join(REPO, 'cache', 'consensus_cache.pkl')
os.makedirs(FIG, exist_ok=True)

COLORS = ['#D7301F', '#08519C', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C']
np.random.seed(42)

print("Loading cached results...")
with open(CACHE, 'rb') as f:
    d = pickle.load(f)
results = d['results']
best_k = d['best_k']
labels = d['labels']
os_df = d['os_df']

n_subtypes = len(set(labels))
consensus = results[best_k]['consensus']
order = np.argsort(labels)

# ---------- panel C and D stats ----------
s1_os = os_df[os_df['subtype'] == 0]
s2_os = os_df[os_df['subtype'] == 1]
lr_tcga = logrank_test(s1_os['os_time'] / 365.25, s2_os['os_time'] / 365.25,
                       s1_os['os'].astype(int), s2_os['os'].astype(int))

tcga_cox = CoxPHFitter()
tcga_cox_df = os_df[['os_time', 'os', 'subtype']].copy()
tcga_cox_df['os'] = tcga_cox_df['os'].astype(int)
tcga_cox_df['os_time'] = tcga_cox_df['os_time'] / 365.25
tcga_cox.fit(tcga_cox_df, 'os_time', 'os', formula='subtype')
tcga_hr = tcga_cox.summary.loc['subtype', 'exp(coef)']
tcga_hr_lo = tcga_cox.summary.loc['subtype', 'exp(coef) lower 95%']
tcga_hr_hi = tcga_cox.summary.loc['subtype', 'exp(coef) upper 95%']

stages = ['I', 'II', 'III', 'IV']
total = len(os_df)
s1_all, s2_all = [], []
for st in stages:
    sub = os_df[os_df['stage'] == st]
    s1_all.append((sub['subtype'] == 0).sum() / total * 100)
    s2_all.append((sub['subtype'] == 1).sum() / total * 100)

print(f"Figure 2 stats: best_k={best_k}, n_subtypes={n_subtypes}, n={total}")
print(f"  logrank P = {lr_tcga.p_value:.3e}")
print(f"  Cox HR = {tcga_hr:.3f} (95% CI {tcga_hr_lo:.3f}-{tcga_hr_hi:.3f})")
print(f"  S1 n = {int((labels == 0).sum())}, "
      f"S2 n = {int((labels == 1).sum())}")

# ---------- figure ----------
FS_PANEL = 10; FS_LABEL = 8; FS_TICK = 7; FS_LEG = 7
FS_CBAR = 7; FS_KM_P = 7; FS_KM_HR = 6.5; FS_BAR = 6; FS_SLAB = 11

fig = plt.figure(figsize=(6.71, 5.01))
fig.patch.set_facecolor('white')

W = 0.40; H_TOP = 0.40; H_BOT = 0.40
LEFT_X = 0.04; RIGHT_X = 0.55

# A: consensus matrix
ax1 = fig.add_axes([LEFT_X, 0.57, W, H_TOP])
ax1.set_facecolor('white')
im = ax1.imshow(consensus[order][:, order], cmap='RdPu',
                vmin=0, vmax=1, interpolation='antialiased', rasterized=True)
for s in range(n_subtypes):
    pos = np.sum(labels[order] <= s)
    ax1.axhline(y=pos - 0.5, color='white', lw=0.5)
    ax1.axvline(x=pos - 0.5, color='white', lw=0.5)
    lo = 0 if s == 0 else int(np.sum(labels[order] <= s - 1))
    hi = int(pos)
    cx = cy = (lo + hi) / 2
    col = '#B71C1C' if s == 0 else '#0D47A1'
    ax1.text(cx, cy, f'S{s+1}', ha='center', va='center', fontsize=FS_SLAB,
             fontweight='bold', color=col,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       alpha=0.85, edgecolor='none'))

ax1.text(-0.04, 1.08, 'A', transform=ax1.transAxes,
         fontsize=FS_PANEL, fontweight='bold', va='top')
for spine in ax1.spines.values():
    spine.set_linewidth(0.8)
cbar = plt.colorbar(im, ax=ax1, shrink=0.8, fraction=0.10, pad=0.04)
cbar.set_label('Consensus index', fontsize=FS_LABEL, fontweight='bold')
cbar.ax.tick_params(labelsize=FS_TICK)
ax1.grid(False)
ax1.tick_params(labelsize=FS_TICK)
ax1.set_xlabel('Samples (ordered by subtype assignment)',
               fontsize=FS_LABEL, fontweight='bold')
ax1.set_ylabel('Samples (ordered by subtype assignment)',
               fontsize=FS_LABEL, fontweight='bold')

# B: metrics
ax2 = fig.add_axes([RIGHT_X, 0.57, W, H_TOP])
ax2.set_facecolor('white')
ks = sorted(results.keys())
sils = [results[k]['silhouette'] for k in ks]
cophs = [results[k]['cophenetic'] for k in ks]
pacs = [results[k]['pac'] for k in ks]
ax2.plot(ks, sils, 'o-', color='#E74C3C', lw=1.2, ms=4,
         label='Silhouette score')
ax2_twin = ax2.twinx()
ax2_twin.plot(ks, cophs, 's-', color='#3498DB', lw=1.2, ms=4,
              label='Cophenetic correlation')
ax2_twin.plot(ks, pacs, '^--', color='#E67E22', lw=1.0, ms=4, label='PAC')
ax2.set_xlabel('Number of clusters (k)',
               fontsize=FS_LABEL, fontweight='bold')
ax2.set_ylabel('Silhouette score',
               fontsize=FS_LABEL, fontweight='bold', color='black')
ax2_twin.set_ylabel('Cophenetic correlation / PAC',
                    fontsize=FS_LABEL, fontweight='bold', color='black')
ax2.text(0, 1.08, 'B', transform=ax2.transAxes,
         fontsize=FS_PANEL, fontweight='bold', va='top')
ax2.set_xticks(ks)
ax2.grid(False)
l1, lab1 = ax2.get_legend_handles_labels()
l2, lab2 = ax2_twin.get_legend_handles_labels()
ax2.legend(l1 + l2, lab1 + lab2, fontsize=FS_LEG, loc='center right',
           bbox_to_anchor=(1.0, 0.7), framealpha=0.9)
for spine in ax2.spines.values():
    spine.set_linewidth(0.8)
ax2.tick_params(labelsize=FS_TICK)
ax2_twin.tick_params(labelsize=FS_TICK)

# C: KM
ax3 = fig.add_axes([LEFT_X, 0.05, W, H_BOT])
kmf = KaplanMeierFitter()
for s in range(n_subtypes):
    sub = os_df[os_df['subtype'] == s]
    if len(sub) > 0:
        kmf.fit(sub['os_time'] / 365.25, sub['os'].astype(int),
                label=f'S{s+1} (n={len(sub)})')
        kmf.plot_survival_function(ax=ax3, color=COLORS[s], lw=1.2)

p_str = (f'P = {lr_tcga.p_value:.1e}' if lr_tcga.p_value < 0.01
         else f'P = {lr_tcga.p_value:.3f}')
ax3.text(0.95, 0.20, p_str, transform=ax3.transAxes, fontsize=FS_KM_P,
         fontweight='bold', ha='right',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                   alpha=0.9, edgecolor='#999'))
hr_str = f'HR = {tcga_hr:.2f} (95% CI: {tcga_hr_lo:.2f}\u2013{tcga_hr_hi:.2f})'
ax3.text(0.95, 0.10, hr_str, transform=ax3.transAxes, fontsize=FS_KM_HR,
         fontweight='bold', ha='right', color='#555')
ax3.set_xlabel('Time (years)', fontsize=FS_LABEL, fontweight='bold')
ax3.set_ylabel('Survival probability',
               fontsize=FS_LABEL, fontweight='bold')
ax3.text(-0.04, 1.08, 'C', transform=ax3.transAxes,
         fontsize=FS_PANEL, fontweight='bold', va='top')
ax3.grid(False)
ax3.tick_params(labelsize=FS_TICK)
ax3.legend(fontsize=FS_LEG)
for spine in ax3.spines.values():
    spine.set_linewidth(0.8)

# D: stage distribution
ax4 = fig.add_axes([RIGHT_X, 0.05, W, H_BOT])
C1_S1, C2_S2 = '#B71C1C', '#0D47A1'
x = np.arange(len(stages))
width = 0.55
ax4.bar(x, s1_all, width, color=C1_S1, alpha=0.85, label='S1')
ax4.bar(x, s2_all, width, bottom=s1_all, color=C2_S2, alpha=0.85, label='S2')
for i in range(len(stages)):
    s1v, s2v = s1_all[i], s2_all[i]
    ax4.text(i, s1v / 2, f'{s1v:.1f}%', ha='center', va='center',
             color='white', fontsize=FS_BAR, fontweight='bold')
    ax4.text(i, s1v + s2v / 2, f'{s2v:.1f}%', ha='center', va='center',
             color='white', fontsize=FS_BAR, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels([f'Stage {st}' for st in stages],
                    fontsize=FS_TICK, fontweight='bold')
ax4.set_xlabel('Tumor stage', fontsize=FS_LABEL, fontweight='bold')
ax4.set_ylabel('Proportion with stage (%)',
               fontsize=FS_LABEL, fontweight='bold')
ax4.text(0, 1.08, 'D', transform=ax4.transAxes,
         fontsize=FS_PANEL, fontweight='bold', va='top')
ax4.set_ylim(0, 65)
ax4.grid(False)
ax4.tick_params(axis='y', labelsize=FS_TICK)
ax4.tick_params(axis='x', labelsize=FS_TICK)
ax4.legend(fontsize=FS_LEG, loc='upper right')
for spine in ax4.spines.values():
    spine.set_linewidth(0.8)

output = os.path.join(FIG, 'Figure_2_consensus_clustering_v7_pubsize.png')
fig.savefig(output, dpi=600, bbox_inches='tight', facecolor='white')
svg_out = os.path.join(FIG, 'Figure_2_consensus_clustering_v7_pubsize.svg')
fig.savefig(svg_out, format='svg', bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {output}")
print(f"Saved: {svg_out}")