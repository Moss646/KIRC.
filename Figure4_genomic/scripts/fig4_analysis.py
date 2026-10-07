# fig4_analysis.py
# Figure 4 (genomic) — analysis step.
# Reads raw TCGA-KIRC MAF + GISTIC CNV, writes intermediate CSVs to data/processed/.

import io
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

import fig4_params as P

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')


def main():
    sub = pd.read_csv(os.path.join(P.PROC_DIR, 'subtype_assignment_balanced.csv'))
    sub_map = dict(zip(sub['Patient'], sub['Subtype']))

    # ---- Panel A: mutation waterfall ----
    maf = pd.read_csv(os.path.join(P.RAW_DIR, 'TCGA-KIRC.merged.maf'),
                      sep='\t', comment='#', low_memory=False)
    maf = maf[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification']].copy()
    maf = maf[maf['Variant_Classification'].isin(P.NONSYN)].copy()
    maf['Patient'] = maf['Tumor_Sample_Barcode'].str[:12]
    print(f'  MAF after nonsynonymous filter: {len(maf)} variants')

    maf_pats = sorted(set(maf['Patient']) & set(sub_map))
    s1_mut_pats = sorted([p for p in maf_pats if sub_map[p] == 'Subtype_1'])
    s2_mut_pats = sorted([p for p in maf_pats if sub_map[p] == 'Subtype_2'])
    n_s1_mut, n_s2_mut = len(s1_mut_pats), len(s2_mut_pats)
    mut_pats_ordered = s1_mut_pats + s2_mut_pats

    rows = []
    onco = []
    for gene in P.MUT_GENES_ORDER:
        s1_n = len(maf[(maf['Hugo_Symbol'] == gene) &
                       (maf['Patient'].isin(s1_mut_pats))]['Patient'].unique())
        s2_n = len(maf[(maf['Hugo_Symbol'] == gene) &
                       (maf['Patient'].isin(s2_mut_pats))]['Patient'].unique())
        s1_no = n_s1_mut - s1_n
        s2_no = n_s2_mut - s2_n
        p = fisher_exact([[s1_n, s1_no], [s2_n, s2_no]])[1] if (s1_n + s2_n) > 0 else 1.0

        rows.append({'gene': gene, 'S1_mut_n': s1_n, 'S2_mut_n': s2_n,
                     'n_s1': n_s1_mut, 'n_s2': n_s2_mut, 'p_raw': p})
        print(f'  Panel A Fisher: {gene} S1={s1_n}/{n_s1_mut} '
              f'S2={s2_n}/{n_s2_mut} p={p:.4e}')

        for pat in mut_pats_ordered:
            g = maf[(maf['Hugo_Symbol'] == gene) & (maf['Patient'] == pat)]
            if len(g) > 0:
                onco.append({'gene': gene, 'patient': pat,
                             'variant_class': g.iloc[0]['Variant_Classification']})

    mut_stats = pd.DataFrame(rows)
    mut_stats['p_bh'] = P.bh_correct(mut_stats['p_raw'].values)
    mut_stats['star'] = mut_stats['p_bh'].apply(P.pval_to_stars)

    pd.DataFrame({
        'patient': mut_pats_ordered,
        'subtype': ['Subtype_1'] * n_s1_mut + ['Subtype_2'] * n_s2_mut,
    }).to_csv(os.path.join(P.PROC_DIR, 'fig4_panelA_samples.csv'), index=False)

    mut_stats.to_csv(os.path.join(P.PROC_DIR, 'fig4_mutation_stats.csv'), index=False)
    pd.DataFrame(onco).to_csv(os.path.join(P.PROC_DIR, 'fig4_oncoprint.csv'), index=False)

    print('  Panel A BH correction (7 tests):')
    for _, r in mut_stats.iterrows():
        print(f"    {r['gene']:8s}  p={r['p_raw']:.4f} -> BH={r['p_bh']:.4f} {r['star']}")

    # ---- Panel B: CNV frequency ----
    cnv = pd.read_csv(os.path.join(P.RAW_DIR, 'TCGA-KIRC.gistic.tsv'),
                      sep='\t', index_col=0)
    cnv.columns = [c[:12] for c in cnv.columns]

    cnv_pats = sorted(set(cnv.columns) & set(sub_map))
    s1_cnv_pats = sorted([p for p in cnv_pats if sub_map[p] == 'Subtype_1'])
    s2_cnv_pats = sorted([p for p in cnv_pats if sub_map[p] == 'Subtype_2'])
    n_s1_cnv, n_s2_cnv = len(s1_cnv_pats), len(s2_cnv_pats)

    cnv_rows = []
    for grp, genes in P.CNV_GROUPS.items():
        for g in genes:
            if g not in cnv.index:
                continue

            s1v = cnv.loc[g, s1_cnv_pats]
            s2v = cnv.loc[g, s2_cnv_pats]

            s1_del, s2_del = int((s1v < 0).sum()), int((s2v < 0).sum())
            s1_amp, s2_amp = int((s1v > 0).sum()), int((s2v > 0).sum())

            p_del = fisher_exact([[s1_del, n_s1_cnv - s1_del],
                                  [s2_del, n_s2_cnv - s2_del]])[1] \
                if (s1_del + s2_del) > 0 else 1.0
            p_amp = fisher_exact([[s1_amp, n_s1_cnv - s1_amp],
                                  [s2_amp, n_s2_cnv - s2_amp]])[1] \
                if (s1_amp + s2_amp) > 0 else 1.0

            cnv_rows.append({
                'gene': g, 'group': grp,
                'S1_del_pct': s1_del / n_s1_cnv * 100,
                'S2_del_pct': s2_del / n_s2_cnv * 100,
                'S1_amp_pct': s1_amp / n_s1_cnv * 100,
                'S2_amp_pct': s2_amp / n_s2_cnv * 100,
                'p_del': p_del, 'p_amp': p_amp,
            })

    cnv_df = pd.DataFrame(cnv_rows)

    all_p = cnv_df['p_del'].tolist() + cnv_df['p_amp'].tolist()
    adj = P.bh_correct(all_p)
    n = len(cnv_df)
    cnv_df['p_del_adj'] = adj[:n]
    cnv_df['p_amp_adj'] = adj[n:]
    cnv_df['star_del'] = cnv_df['p_del_adj'].apply(P.pval_to_stars)
    cnv_df['star_amp'] = cnv_df['p_amp_adj'].apply(P.pval_to_stars)

    pd.DataFrame({
        'patient': s1_cnv_pats + s2_cnv_pats,
        'subtype': ['Subtype_1'] * n_s1_cnv + ['Subtype_2'] * n_s2_cnv,
    }).to_csv(os.path.join(P.PROC_DIR, 'fig4_panelB_samples.csv'), index=False)

    cnv_df.to_csv(os.path.join(P.PROC_DIR, 'fig4_cnv_stats.csv'), index=False)

    print('  BH correction applied (28 tests: 14 genes x 2 directions)')
    for _, r in cnv_df.iterrows():
        print(f"    {r['gene']:10s}  del p={r['p_del']:.4f}->adj={r['p_del_adj']:.4f}"
              f"  amp p={r['p_amp']:.4f}->adj={r['p_amp_adj']:.4f}")

    print('\nAnalysis complete. Intermediate CSVs written to data/processed/')


if __name__ == '__main__':
    import traceback
    exit_code = 0
    try:
        main()
    except SystemExit as e:
        exit_code = e.code if isinstance(e.code, int) else 1
    except BaseException:
        traceback.print_exc()
        with open("error.log", "w", encoding="utf-8") as f:
            traceback.print_exc(file=f)
        print("\n[致命错误] 详细信息已写入当前目录下的 error.log")
        exit_code = 1
    finally:
        print("\n" + "=" * 60)
        print("  运行结束（或失败）。按任意键退出...")
        print("=" * 60)
        if os.name == "nt":
            import msvcrt
            msvcrt.getch()
        else:
            input()
    sys.exit(exit_code)