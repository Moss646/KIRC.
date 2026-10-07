#!/usr/bin/env python3
# build_classifier.py
# Train the NCC subtype classifier on TCGA-KIRC discovery cohort.
# Saves classifier JSON, self-validation figure, and ICGC validation skeleton.

import sys, os, io, atexit, time, warnings, pickle, json, re
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
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
DATA = os.path.join(REPO, 'data')
FIG = os.path.join(REPO, 'results')
BG = '#FAFAFA'
os.makedirs(FIG, exist_ok=True)


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
    m = re.match(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})', s)
    return m.group(1) if m else s[:12]


def _norm_os(x):
    s = str(x).strip().lower()
    if s in ('1', 'true', 'dead', 'deceased'):
        return 1
    if s in ('0', 'false', 'alive', 'living'):
        return 0
    try:
        return int(float(s))
    except Exception:
        return None


def _load_clinical():
    info_path = os.path.join(DATA, 'TCGA-KIRC_clinical_info.csv')
    surv_path = os.path.join(DATA, 'TCGA-KIRC_clinical_survival.csv')
    for p in (info_path, surv_path):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    info = pd.read_csv(info_path)
    surv = pd.read_csv(surv_path)

    sid = _pick(surv, ['bcr_patient_barcode', 'submitter_id',
                       'case_id', 'sample'])
    sos = _pick(surv, ['OS', 'vital_status'])
    stime = _pick(surv, ['OS.time', 'os_time_days',
                         'days_to_death', 'days_to_last_follow_up'])

    surv2 = pd.DataFrame({
        'Patient':      surv[sid].map(_to_patient),
        'os':           surv[sos].map(_norm_os),
        'os_time_days': pd.to_numeric(surv[stime], errors='coerce'),
    }).drop_duplicates('Patient').set_index('Patient')
    return surv2.dropna(subset=['os', 'os_time_days'])


# ---------- load data ----------
EXPR = pd.read_csv(os.path.join(DATA, 'KIRC_expr_log2_tpm.csv'), index_col=0)
CLIN = _load_clinical()

with open(os.path.join(DATA, 'balanced_genes.txt')) as f:
    sel_genes = [l.strip() for l in f if l.strip()]
sel_genes = [g for g in sel_genes if g in EXPR.index]
print(f"Genes: {len(sel_genes)}, Patients: {len(EXPR.columns)}")

# ---------- build classifier from full TCGA ----------
X = EXPR.loc[sel_genes].T.values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

sub = pd.read_csv(os.path.join(DATA, 'subtype_assignment_balanced.csv'))
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
labels = np.array([
    1 if sub_map.get(s, 'Subtype_2') == 'Subtype_1' else 0
    for s in EXPR.columns
])

centroid_s1 = X_scaled[labels == 1].mean(axis=0)
centroid_s2 = X_scaled[labels == 0].mean(axis=0)

classifier = {
    'genes': sel_genes,
    'n_genes': len(sel_genes),
    'centroid_s1': centroid_s1.tolist(),
    'centroid_s2': centroid_s2.tolist(),
    'scaler_mean': scaler.mean_.tolist(),
    'scaler_std': scaler.scale_.tolist() if hasattr(scaler, 'scale_')
                  else scaler.var_.tolist(),
    'subtype_names': {0: 'S2', 1: 'S1'},
    'tcga_performance': {
        'S1': int((labels == 1).sum()),
        'S2': int((labels == 0).sum()),
    }
}

with open(os.path.join(DATA, 'subtype_classifier.json'), 'w') as f:
    json.dump(classifier, f, indent=2)
print(f"Classifier saved: {os.path.join(DATA, 'subtype_classifier.json')}")

# ---------- self-validation ----------
dist_s1 = np.linalg.norm(X_scaled - centroid_s1, axis=1)
dist_s2 = np.linalg.norm(X_scaled - centroid_s2, axis=1)
pred_labels = (dist_s1 < dist_s2).astype(int)
acc = accuracy_score(labels, pred_labels)
print(f"\nSelf-validation accuracy: {acc:.1%}")

discordant = labels != pred_labels
print(f"Discordant patients: {discordant.sum()}/{len(labels)} "
      f"({discordant.sum()/len(labels)*100:.1f}%)")

confidence = np.abs(dist_s1 - dist_s2)
high_conf = confidence > np.median(confidence)
low_conf = ~high_conf

# ---------- ICGC R skeleton ----------
r_code = '''# ICGC RECA-EU external validation
# 1. download donor/specimen/exp_seq from
#    https://dcc.icgc.org/releases/current/Projects/RECA-EU
# 2. put this script and subtype_classifier.json in the same dir
# 3. run in R

library(jsonlite)
library(survival)
library(survminer)

classifier <- fromJSON("subtype_classifier.json")

donor <- read.delim("donor.RECA-EU.tsv.gz")
specimen <- read.delim("specimen.RECA-EU.tsv.gz")
expr <- read.delim("exp_seq.RECA-EU.tsv.gz")

tumor_specimens <- specimen[specimen$specimen_type == "Primary tumour - solid tissue", ]

library(reshape2)
expr_tumor <- expr[expr$icgc_specimen_id %in% tumor_specimens$icgc_specimen_id, ]
expr_mat <- dcast(expr_tumor, gene_id ~ icgc_specimen_id,
                  value.var = "normalized_read_count", fun.aggregate = mean)
rownames(expr_mat) <- expr_mat$gene_id
expr_mat$gene_id <- NULL

library(org.Hs.eg.db)
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(expr_mat),
                       column = "SYMBOL", keytype = "ENSEMBL")
expr_mat$symbol <- gene_symbols
expr_mat <- expr_mat[!is.na(expr_mat$symbol), ]
expr_by_symbol <- aggregate(. ~ symbol, data = expr_mat, FUN = mean)
rownames(expr_by_symbol) <- expr_by_symbol$symbol
expr_by_symbol$symbol <- NULL

common_genes <- intersect(classifier$genes, rownames(expr_by_symbol))
expr_sub <- as.matrix(expr_by_symbol[common_genes, ])
expr_log <- log2(expr_sub + 1)

expr_z <- scale(t(expr_log), center = classifier$scaler_mean,
                scale = classifier$scaler_std)

dist_s1 <- sqrt(rowSums((expr_z - classifier$centroid_s1)^2))
dist_s2 <- sqrt(rowSums((expr_z - classifier$centroid_s2)^2))
subtypes <- ifelse(dist_s1 < dist_s2, "S1", "S2")

cat("ICGC external validation completed!\\n")
cat(sprintf("Predicted subtypes: S1=%d, S2=%d\\n",
            sum(subtypes == "S1"), sum(subtypes == "S2")))
'''

with open(os.path.join(DATA, 'icgc_validation.R'), 'w') as f:
    f.write(r_code)
print(f"ICGC R code: {os.path.join(DATA, 'icgc_validation.R')}")

# ---------- figure ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.patch.set_facecolor(BG)

# A: confidence distribution
ax = axes[0]; ax.set_facecolor(BG)
for label, color, name in [(1, '#B71C1C', 'S1'), (0, '#0D47A1', 'S2')]:
    conf_sub = confidence[labels == label]
    ax.hist(conf_sub, bins=30, color=color, alpha=0.6, label=name,
            edgecolor='white')
ax.axvline(x=np.median(confidence), color='#333', ls='--', lw=1.5)
ax.set_xlabel('Classification Confidence', fontsize=10)
ax.set_ylabel('Patients', fontsize=10)
ax.set_title('A. Classification Confidence', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15, axis='y')

# B: KM by confidence
ax = axes[1]; ax.set_facecolor(BG)
kmf = KaplanMeierFitter()
for mask, label, color in [(high_conf, 'High Confidence', '#2ECC71'),
                            (low_conf, 'Low Confidence', '#E74C3C')]:
    idx = np.where(mask)[0]
    rows = []
    for i in idx:
        pid = EXPR.columns[i]
        if pid in CLIN.index:
            r = CLIN.loc[pid]
            rows.append({'time': r['os_time_days'],
                         'status': r['os']})
    df = pd.DataFrame(rows).dropna()
    if len(df) > 0:
        kmf.fit(df['time'] / 365.25, df['status'].astype(int),
                label=f'{label} (n={len(df)})')
        kmf.plot_survival_function(ax=ax, color=color, lw=2)

# log-rank between high/low confidence
high_rows = [{'time': CLIN.loc[EXPR.columns[i], 'os_time_days'],
              'status': CLIN.loc[EXPR.columns[i], 'os']}
             for i in np.where(high_conf)[0] if EXPR.columns[i] in CLIN.index]
low_rows = [{'time': CLIN.loc[EXPR.columns[i], 'os_time_days'],
             'status': CLIN.loc[EXPR.columns[i], 'os']}
            for i in np.where(low_conf)[0] if EXPR.columns[i] in CLIN.index]
df_h = pd.DataFrame(high_rows).dropna()
df_l = pd.DataFrame(low_rows).dropna()
lr_conf = logrank_test(df_h['time'], df_l['time'],
                       df_h['status'].astype(int),
                       df_l['status'].astype(int))

ax.set_xlabel('Time (years)', fontsize=10)
ax.set_ylabel('Overall Survival', fontsize=10)
ax.set_title(f'B. Survival by Confidence\nLog-rank p={lr_conf.p_value:.2e}',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.15)
ax.legend(fontsize=8)

# C: centroids
ax = axes[2]; ax.set_facecolor(BG)
centroids = np.vstack([centroid_s2, centroid_s1])
im = ax.imshow(centroids, aspect='auto', cmap='RdBu_r', vmin=-1.5, vmax=1.5)
ax.set_yticks([0, 1])
ax.set_yticklabels(['S2', 'S1'], fontsize=10)
ax.set_xticks([])
top_idx = np.argsort(np.abs(centroid_s1 - centroid_s2))[-10:]
for idx in top_idx:
    ax.annotate(sel_genes[idx], (idx, -1.5), fontsize=5,
                rotation=90, ha='right', va='top')
ax.set_title('C. Subtype Centroids (top genes)',
             fontsize=11, fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.6, label='Z-score')

fig.suptitle('Subtype Classifier: Construction & Validation',
             fontsize=13, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIG, 'classifier_validation.png'),
            dpi=200, bbox_inches='tight', facecolor=BG)
plt.close()
print(f"classifier_validation.png -> "
      f"{os.path.join(FIG, 'classifier_validation.png')}")
print(f"\nClassifier ready for external validation!")
print(f"  {os.path.join(DATA, 'subtype_classifier.json')}")
print(f"  {os.path.join(DATA, 'icgc_validation.R')}")