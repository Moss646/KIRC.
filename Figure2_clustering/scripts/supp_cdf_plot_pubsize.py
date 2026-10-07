# supp_cdf_plot_pubsize.py
# Publication-size CDF figure (180 x 90 mm, 600 dpi).

import sys, os, io, atexit, pickle, warnings
import numpy as np
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
plt.rcParams['svg.fonttype'] = 'none'

WORK = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(WORK)
CACHE = os.path.join(REPO, 'cache', 'consensus_cache.pkl')
FIG = os.path.join(REPO, 'results')
os.makedirs(FIG, exist_ok=True)

with open(CACHE, 'rb') as f:
    d = pickle.load(f)
results = d['results']

cdf_data = {}
for k in sorted(results.keys()):
    consensus = results[k]['consensus']
    upper = consensus[np.triu_indices_from(consensus, k=1)]
    upper = upper[upper > 0]
    if len(upper) > 0:
        cdf_x = np.sort(upper)
        cdf_y = np.arange(1, len(cdf_x) + 1) / len(cdf_x)
        cdf_auc = np.trapezoid(cdf_y, cdf_x)
    else:
        cdf_x, cdf_y = np.array([0, 1]), np.array([0, 1])
        cdf_auc = 0.5
    cdf_data[k] = {'cdf_x': cdf_x, 'cdf_y': cdf_y, 'cdf_auc': cdf_auc}

ks = sorted(cdf_data.keys())
colors = ['#D7301F', '#08519C', '#2ECC71', '#F39C12', '#9B59B6']

FS_PANEL = 10; FS_LABEL = 8; FS_TICK = 7; FS_LEG = 7; FS_BAR = 7

fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.54))
fig.patch.set_facecolor('white')

# A: CDF curves
ax = axes[0]
ax.set_facecolor('white')
for i, k in enumerate(ks):
    ax.plot(cdf_data[k]['cdf_x'], cdf_data[k]['cdf_y'],
            color=colors[i], lw=1.2, label=f'k = {k}', alpha=0.85)
ax.set_xlabel('Consensus index', fontsize=FS_LABEL, fontweight='bold')
ax.set_ylabel('CDF', fontsize=FS_LABEL, fontweight='bold')
ax.text(-0.06, 1.08, 'A', transform=ax.transAxes,
        fontsize=FS_PANEL, fontweight='bold', va='top')
ax.tick_params(labelsize=FS_TICK)
ax.legend(fontsize=FS_LEG, loc='upper left', framealpha=0.9)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
for spine in ax.spines.values():
    spine.set_linewidth(0.8)

# B: delta AUC
ax2 = axes[1]
ax2.set_facecolor('white')
aucs = [cdf_data[k]['cdf_auc'] for k in ks]
deltas = [aucs[i] - aucs[i - 1] for i in range(1, len(aucs))]
delta_ks = ks[1:]
bars = ax2.bar(delta_ks, deltas, color='#87CEEB', alpha=0.85,
               edgecolor='#5DADE2', linewidth=0.8, width=0.6)
ax2.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
ax2.set_xticks(delta_ks)
ax2.set_xticklabels([f'k = {k}' for k in delta_ks],
                    fontsize=FS_TICK, fontweight='bold')
ax2.set_xlabel('Number of clusters (k)',
               fontsize=FS_LABEL, fontweight='bold')
ax2.set_ylabel('Change in AUC (\u0394 AUC)',
               fontsize=FS_LABEL, fontweight='bold')
ax2.text(-0.06, 1.08, 'B', transform=ax2.transAxes,
         fontsize=FS_PANEL, fontweight='bold', va='top')
ax2.tick_params(labelsize=FS_TICK)
for spine in ax2.spines.values():
    spine.set_linewidth(0.8)
for bar, val in zip(bars, deltas):
    ax2.text(bar.get_x() + bar.get_width() / 2., val + 0.002,
             f'{val:.3f}', ha='center', va='bottom',
             fontsize=FS_BAR, fontweight='bold', color='#2C3E50')

plt.tight_layout(w_pad=3)

out_png = os.path.join(FIG, 'Supplementary_Figure_CDF_analysis_pubsize.png')
out_svg = os.path.join(FIG, 'Supplementary_Figure_CDF_analysis_pubsize.svg')
fig.savefig(out_png, dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(out_svg, format='svg', bbox_inches='tight', facecolor='white')
plt.close()

print(f'Saved: {out_png}')
print(f'Saved: {out_svg}')
print('AUC:', {k: round(cdf_data[k]['cdf_auc'], 4) for k in ks})
print('Delta AUC:', {str(delta_ks[i]): round(deltas[i], 4)
                     for i in range(len(deltas))})