# Figure6_associations_lipid_immune

Lipid-metabolism x immune-feature correlation analysis. 30 genes
(14 fatty-acid metabolism + 16 phospholipid-remodeling) x 7 immune
features = 210 pairs. Spearman correlation across four analytic schemes:
full cohort, partial correlation adjusted for immune + stromal content,
and within-subtype correlation.

## Inputs

| File | Source | How to get |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/tcga_kirc/` |
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/` |
| `TIDE_official_results.csv` | Figure7_therapeutic_analysis_pubsize | copy to `data/` |
| `xcell_scores_raw.csv` | Figure5_immune | copy to `data/` |

Note: `KIRC_expr_log2_tpm.csv` lives under `data/tcga_kirc/`; the other
three files live directly under `data/`. When copying `xcell_scores_raw.csv`
from `Figure5_immune/data/xcell/`, place it directly under
`Figure6_associations_lipid_immune/data/` without the `xcell/` subdirectory.

## Pipeline steps

1. `metab_immune_corr_v4.py` - merge expression, subtypes, TIDE, xCell,
   compute Spearman, partial Spearman (adjusted for ImmuneScore, then
   ImmuneScore + StromaScore), and subtype-stratified Spearman for all
   210 gene x feature pairs. BH correction within each scheme.
2. `fig_Figure6_main_v1.py` - render the main figure (Panel A: heatmap of
   adjusted rho; Panel B: core association dot plot across four strata).

## Outputs

| File | Location | Consumer |
|---|---|---|
| `metabolic_immune_corr_30genes_v1.csv` | `results/` | internal |
| `Figure6_metab_immune_corr_v1.{png,svg}` | `results/` | main text |

## Path configuration

`config.py` resolves the expression matrix as
`<MEL_DATA_ROOT>/tcga_kirc/KIRC_expr_log2_tpm.csv`, and the other three
inputs as `<PROJECT>/data/`. The provided `run_all.py` sets
`MEL_DATA_ROOT` to `<PROJECT>/data/` automatically.

## Run

    cd Figure6_associations_lipid_immune/scripts
    python run_all.py

Runtime: under 1 minute.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).