# Figure4_genomic

Mutation landscape (driver MAF) and copy-number alteration (GISTIC) of
the two subtypes, plus 3p driver deletion rates and FAO9 score by
mutation status.

## Inputs

| File | Source | How to get |
|---|---|---|
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/processed/` |
| `TCGA-KIRC.merged.maf` | GDC | download to `data/raw/` |
| `TCGA-KIRC.gistic.tsv` | GDAC Firehose | download to `data/raw/` |
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/raw/` |
| `nonncc9_fao_score_per_sample_v1.csv` | Figure8_External_validation_OS_KM_score_combined | copy to `data/raw/` |

The MAF and GISTIC files come from the TCGA-KIRC GDAC Firehose legacy
release; the FAO9 score table is produced by Figure8's
`nonncc9_fao_cross_cohort_v1.py`.

## Pipeline steps

1. `make_fao9_input.py` - convert the Figure8 FAO9 score table into the
   two-column TCGA-only table expected by the driver-mutation scripts.
2. `fig4_analysis.py` - parse MAF and GISTIC, compute Fisher exact tests
   per driver gene, apply BH correction, write intermediate CSVs.
3. `cnv_lipid_remodeling_16genes_check_v1.py` - CNA check for the 16
   phospholipid remodeling genes (BH pooled across the 32 deletion/gain
   tests).
4. `driver_mut_fao_v1.py` - driver mutations vs subtype and FAO9 score,
   plus 3p driver deletion rates.
5. `fig_driver_3p_del_FAO9_v1.py` - two-panel supplementary figure.
6. `Fig4_plot_v30_pubsize.py` - assemble the main Figure 4.

## Outputs

| File | Location | Consumer |
|---|---|---|
| `fao9_immune_input_tcga_v1.csv` | `data/processed/` | internal |
| `fig4_panelA_samples.csv` | `data/processed/` | internal |
| `fig4_mutation_stats.csv` | `data/processed/` | internal |
| `fig4_oncoprint.csv` | `data/processed/` | internal |
| `fig4_panelB_samples.csv` | `data/processed/` | internal |
| `fig4_cnv_stats.csv` | `data/processed/` | internal |
| `driver_mut_fao_v1.csv` | `results/tables/` | internal |
| `driver_gistic_del_v1.csv` | `results/tables/` | internal |
| `cnv_lipid_remodeling_16genes_check_v2.csv` | `results/tables/` | supplement |
| `Figure_4_genomic_v31_pubsize.{png,svg}` | `results/figures/` | main text |
| `Supplementary_Fig_driver_3p_del_FAO9_v1.{png,svg}` | `results/figures/` | supplement |

## Run

    cd Figure4_genomic/scripts
    python run_all.py

Note: three scripts (`fig4_analysis.py`, `driver_mut_fao_v1.py`,
`fig_driver_3p_del_FAO9_v1.py`) prompt for a key press on exit. The
`run_all.py` wrapper sets `NOPAUSE=1` for the first one, but if you run
them individually you will need to press a key.

Runtime: 2-5 min.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).