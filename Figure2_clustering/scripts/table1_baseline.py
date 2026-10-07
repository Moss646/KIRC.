# table1_baseline.py
# Table 1: baseline characteristics by NCC subtype (TCGA-KIRC).

import sys, os, io, atexit, re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, chi2_contingency


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

WORK = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(WORK)
DATA = os.path.join(REPO, 'data')
RESULTS = os.path.join(REPO, 'results')
os.makedirs(RESULTS, exist_ok=True)


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


def _norm_stage(x):
    s = str(x).strip().upper()
    m = re.search(r'\b(IV|III|II|I)\b', s)
    return m.group(1) if m else None


def _norm_gender(x):
    s = str(x).strip().lower()
    if s.startswith('m'):
        return 'Male'
    if s.startswith('f'):
        return 'Female'
    return s


def _norm_grade(x):
    s = str(x).strip()
    s = re.sub(r'^G', '', s)
    return s if s else None


info_path = os.path.join(DATA, 'TCGA-KIRC_clinical_info.csv')
if not os.path.exists(info_path):
    raise FileNotFoundError(info_path)

info = pd.read_csv(info_path)

iid = _pick(info, ['bcr_patient_barcode', 'submitter_id',
                   'case_id', 'patient'])
iage = _pick(info, ['age_at_index', 'age'])
igender = _pick(info, ['gender'])
istage = _pick(info, ['ajcc_pathologic_tumor_stage',
                      'ajcc_pathologic_stage', 'tumor_stage', 'stage'])
igrade = _pick(info, ['tumor_grade', 'neoplasm_histologic_grade', 'grade'])

missing = [n for n, v in [
    ('id', iid), ('age', iage), ('gender', igender),
    ('stage', istage), ('grade', igrade),
] if v is None]

if missing:
    print('\n[column match failed] missing:', missing)
    print('info columns:', list(info.columns))
    raise SystemExit(1)

clin = pd.DataFrame({
    'Patient': info[iid].map(_to_patient),
    'age':     pd.to_numeric(info[iage], errors='coerce'),
    'gender':  info[igender].map(_norm_gender),
    'stage':   info[istage].map(_norm_stage),
    'grade':   info[igrade].map(_norm_grade),
}).drop_duplicates('Patient').set_index('Patient')

sub = pd.read_csv(os.path.join(DATA, 'subtype_assignment_balanced.csv'))
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
clin['Subtype'] = clin.index.map(sub_map)
clin = clin.dropna(subset=['Subtype'])

print(f'Table 1 cohort: {len(clin)} patients '
      f'(S1={sum(clin["Subtype"]=="Subtype_1")}, '
      f'S2={sum(clin["Subtype"]=="Subtype_2")})')

rows = [{'Variable': 'Variable',
         'All (n=533)': 'All (n=533)',
         'S1 (n=172)': 'S1 (n=172)',
         'S2 (n=361)': 'S2 (n=361)',
         'P value': 'P value'}]

s1_age = clin[clin['Subtype'] == 'Subtype_1']['age'].dropna()
s2_age = clin[clin['Subtype'] == 'Subtype_2']['age'].dropna()
_, age_p = mannwhitneyu(s1_age, s2_age)
all_age = clin['age'].dropna()

rows.append({
    'Variable': 'Age (years), median (IQR)',
    'All (n=533)': f'{all_age.median():.0f} '
                   f'({all_age.quantile(0.25):.0f}\u2013{all_age.quantile(0.75):.0f})',
    'S1 (n=172)': f'{s1_age.median():.0f} '
                  f'({s1_age.quantile(0.25):.0f}\u2013{s1_age.quantile(0.75):.0f})',
    'S2 (n=361)': f'{s2_age.median():.0f} '
                  f'({s2_age.quantile(0.25):.0f}\u2013{s2_age.quantile(0.75):.0f})',
    'P value': f'{age_p:.4f}',
})

for var, label in [('gender', 'Gender'),
                   ('stage', 'Pathological stage'),
                   ('grade', 'Nuclear grade')]:
    data = clin[['Subtype', var]].dropna()
    s1_dist = data[data['Subtype'] == 'Subtype_1'][var].value_counts()
    s2_dist = data[data['Subtype'] == 'Subtype_2'][var].value_counts()
    all_dist = data[var].value_counts()
    cats = sorted(set(s1_dist.index) | set(s2_dist.index), key=str)

    s1_t = len(data[data['Subtype'] == 'Subtype_1'])
    s2_t = len(data[data['Subtype'] == 'Subtype_2'])
    all_t = len(data)

    for c in cats:
        rows.append({
            'Variable': f'  {c}',
            'All (n=533)': f'{int(all_dist.get(c, 0))} '
                           f'({int(all_dist.get(c, 0)) / all_t * 100:.1f}%)',
            'S1 (n=172)': f'{int(s1_dist.get(c, 0))} '
                          f'({int(s1_dist.get(c, 0)) / s1_t * 100:.1f}%)',
            'S2 (n=361)': f'{int(s2_dist.get(c, 0))} '
                          f'({int(s2_dist.get(c, 0)) / s2_t * 100:.1f}%)',
            'P value': '',
        })

    ct = np.array([[int(s1_dist.get(c, 0)), int(s2_dist.get(c, 0))]
                   for c in cats])
    _, chi_p, _, _ = chi2_contingency(ct)
    rows.append({'Variable': label,
                 'All (n=533)': '', 'S1 (n=172)': '',
                 'S2 (n=361)': '', 'P value': f'{chi_p:.2e}'})

df = pd.DataFrame(rows)
out = os.path.join(RESULTS, 'Table1_baseline_characteristics.csv')
df.to_csv(out, index=False)
print(f'Saved: {out}')
print(df.to_string(index=False))