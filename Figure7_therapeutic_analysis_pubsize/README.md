# Figure7_therapeutic_analysis_pubsize

Drug sensitivity prediction (oncoPredict::calcPhenotype with GDSC2
training) and TIDE immune-evasion analysis across the two subtypes.
Panels: A drug sensitivity | B Cohen's d | C residual Cohen's d |
D responder rate | E TIDE scores | F immune markers.

## Inputs

| File | Source | How to get |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/` |
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/` |
| `GDSC2_Expr.rds` | GDSC2 | download to `data/raw/` |
| `GDSC2_Res.rds` | GDSC2 | download to `data/raw/` |

GDSC2 RDS files are available from the GDSC2 download page
(https://www.cancerrxgene.org/downloads/bulk_download).

## Pipeline steps

1. `prepare_oncopredict_inputs.py` - prepare `kirc_expr_matched.csv` and
   `kirc_subtypes.csv` from the expression matrix and subtype assignment.
2. `run_oncoPredict_calcPheno.R` - train GDSC2 ridge/glmnet models and
   predict per-sample IC50 for 9 target drugs.
3. `run_tidepy.py` - run TIDE (tidepy) on the same expression matrix to
   produce immune-evasion scores and responder labels.
4. `Fig7_main_pubsize_v2.py` - render the 6-panel main figure.
5. `Fig7_PanelA_violin_faceted.py` - render the 3x3 faceted violin plot
   of predicted ln(IC50) per drug.

## Outputs

| File | Location | Consumer |
|---|---|---|
| `kirc_expr_matched.csv` | `data/raw/` | internal |
| `kirc_subtypes.csv` | `data/raw/` | internal |
| `oncoPredict_calcPheno_results.csv` | `data/processed/` | internal |
| `oncoPredict_calcPheno_ic50.csv` | `data/processed/` | internal |
| `oncoPredict_calcPheno_<Drug>_ic50.csv` | `data/processed/` | internal |
| `TIDE_official_results.csv` | `data/processed/` | Figure6 |
| `Figure_7_therapeutic_analysis_pubsize_v2.{png,svg}` | `results/figures/` | main text |
| `Fig7_PanelA_violin_faceted.{png,svg}` | `results/figures/` | supplement |

## Run

    cd Figure7_therapeutic_analysis_pubsize/scripts
    python run_all.py

Runtime: 30-120 min. The R oncoPredict step (train + predict for 9
drugs) dominates.

## Notes

- Figure7's main figure does not depend on any external glmnet result
  file. The Panel A drug order is derived downstream from the residual
  Cohen's d, computed inside `Fig7_main_pubsize_v2.py`.
- `run_oncoPredict_calcPheno.R` writes a temporary
  `calcPhenotype_Output/` directory inside `data/processed/`. This is
  cleaned up automatically at the end of the run.
- `tidepy` must be installed in the Python environment: `pip install tidepy`.

## After the run

    # to Figure6
    cp data/processed/TIDE_official_results.csv  ../Figure6_associations_lipid_immune/data/

On Windows, use `copy` instead of `cp`.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).