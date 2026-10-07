# driver_mut_fao_v1.py
# Driver mutations vs subtype and FAO9 score, plus 3p deletion rates.

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig4_params as P

GENES = P.MUT_GENES_ORDER
NONSYN = P.NONSYN


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


sub = pd.read_csv(os.path.join(P.PROC_DIR, 'subtype_assignment_balanced.csv'))
sub['S'] = sub['Subtype'].map({'Subtype_1': 'S1', 'Subtype_2': 'S2'})
sub_map = dict(zip(sub['Patient'], sub['S']))
print(f'subtype map: S1={sum(v=="S1" for v in sub_map.values())}, '
      f'S2={sum(v=="S2" for v in sub_map.values())}')

maf = pd.read_csv(os.path.join(P.RAW_DIR, 'TCGA-KIRC.merged.maf'),
                  sep='\t', comment='#', low_memory=False)
maf = maf[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification']].copy()
maf = maf[maf['Variant_Classification'].isin(NONSYN)].copy()
maf['Patient'] = maf['Tumor_Sample_Barcode'].str[:12]

mut_pats = sorted(set(maf['Patient']) & set(sub_map))
s1_pats = [p for p in mut_pats if sub_map[p] == 'S1']
s2_pats = [p for p in mut_pats if sub_map[p] == 'S2']
print(f'MAF nonsyn variants: {len(maf)}; mutation cohort: '
      f'{len(mut_pats)} (S1={len(s1_pats)}, S2={len(s2_pats)})')

fao9 = pd.read_csv(os.path.join(P.PROC_DIR, 'fao9_immune_input_tcga_v1.csv'))
fao9 = fao9[['sample', 'FAO9']].copy()
fao9['Patient'] = fao9['sample'].str[:12]
fao9_map = dict(zip(fao9['Patient'], fao9['FAO9']))
print(f'FAO9 scores available: {len(fao9_map)}')

rows = []
for gene in GENES:
    g = maf[maf['Hugo_Symbol'] == gene]
    s1_mut = set(g['Patient']) & set(s1_pats)
    s2_mut = set(g['Patient']) & set(s2_pats)
    n1, n2 = len(s1_pats), len(s2_pats)
    m1, m2 = len(s1_mut), len(s2_mut)
    orr, p_fish = fisher_exact([[m1, n1 - m1], [m2, n2 - m2]])

    wt_pats = [p for p in mut_pats if p not in set(g['Patient'])]
    v_mut = [fao9_map[p] for p in mut_pats if p in set(g['Patient']) and p in fao9_map]
    v_wt = [fao9_map[p] for p in wt_pats if p in fao9_map]

    if len(v_mut) >= 3 and len(v_wt) >= 3:
        u, p_mwu = mannwhitneyu(v_mut, v_wt, alternative='two-sided')
        med_mut = pd.Series(v_mut).median()
        med_wt = pd.Series(v_wt).median()
    else:
        p_mwu, med_mut, med_wt = float('nan'), float('nan'), float('nan')

    rows.append(dict(
        gene=gene,
        S1_mut=m1, S1_n=n1, S1_pct=100 * m1 / n1,
        S2_mut=m2, S2_n=n2, S2_pct=100 * m2 / n2,
        OR_S1vsS2=orr, p_fisher=p_fish,
        FAO9_n_mut=len(v_mut), FAO9_n_wt=len(v_wt),
        FAO9_med_mut=med_mut, FAO9_med_wt=med_wt,
        p_mwu_FAO9=p_mwu
    ))

res = pd.DataFrame(rows)
res['p_bh_fisher'] = bh_fdr(res['p_fisher'].tolist())
res['p_bh_mwu'] = bh_fdr(res['p_mwu_FAO9'].tolist())

out_ab = os.path.join(P.PROJ_ROOT, 'results', 'tables', 'driver_mut_fao_v1.csv')
res.to_csv(out_ab, index=False)
print('\n=== A+B: driver mutations vs S1/S2 + FAO9 score ===')
print(res.round(4).to_string(index=False))

gist = pd.read_csv(os.path.join(P.RAW_DIR, 'TCGA-KIRC.gistic.tsv'),
                   sep='\t', index_col=0)

rows_c = []
for gene in ['VHL', 'PBRM1', 'BAP1', 'SETD2']:
    if gene not in gist.index:
        print(f'  [warn] {gene} not in GISTIC table')
        continue
    row = gist.loc[gene]
    samples = gist.columns
    pat = pd.Series(samples, index=samples).str[:12]
    df = pd.DataFrame({'Patient': pat.values, 'val': row.values})
    df = df.groupby('Patient')['val'].min().rename('val').reset_index()
    df = df[df['Patient'].isin(sub_map)].copy()
    df['S'] = df['Patient'].map(sub_map)

    n_all = len(df)
    del_all = int((df['val'] <= -1).sum())
    s1 = df[df['S'] == 'S1']
    s2 = df[df['S'] == 'S2']
    d1, d2 = int((s1['val'] <= -1).sum()), int((s2['val'] <= -1).sum())
    _, p_del = fisher_exact([[d1, len(s1) - d1], [d2, len(s2) - d2]])

    rows_c.append(dict(
        gene=gene, n=n_all,
        del_n=del_all, del_pct=100 * del_all / n_all,
        S1_del=d1, S1_n=len(s1), S1_del_pct=100 * d1 / len(s1),
        S2_del=d2, S2_n=len(s2), S2_del_pct=100 * d2 / len(s2),
        p_fisher_del_S1vsS2=p_del
    ))

res_c = pd.DataFrame(rows_c)
out_c = os.path.join(P.PROJ_ROOT, 'results', 'tables', 'driver_gistic_del_v1.csv')
res_c.to_csv(out_c, index=False)
print('\n=== C: GISTIC gene-level deletion rates (3p drivers) ===')
print(res_c.round(4).to_string(index=False))
print(f'\nSaved: {out_ab}')
print(f'Saved: {out_c}')
import os
print("\n" + "=" * 60)
print("  运行结束（或失败）。按任意键退出...")
print("=" * 60)
if os.name == "nt":
    import msvcrt
    msvcrt.getch()
else:
    input()