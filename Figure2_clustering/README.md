# Figure2_clustering

Consensus clustering (KMeans + average linkage) to define two subtypes,
followed by nearest-centroid classifier construction, NMF validation, and
a survival-blind reproducibility analysis.

## Inputs

| File | Source | How to get |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/` |
| `KIRC_clinical_cleaned.csv` | TCGA-KIRC_prep | copy to `data/` |
| `balanced_genes.txt` | Figure1_gene_selection | copy to `data/` |
| `msigdb_lipid_gene_pool_dedup.csv` | Figure1_gene_selection | copy to `data/` |
| `TCGA-KIRC_clinical_info.csv` | GDC | download to `data/` |
| `TCGA-KIRC_clinical_survival.csv` | GDC | download to `data/` |

The two GDC files are only used by `table1_baseline.py` and by
`compute_clustering.py` for stage / OS linkage. `KIRC_clinical_cleaned.csv`
provides the same data in a cleaner form and is preferred when available.

## Outputs

| File | Location | Consumer |
|---|---|---|
| `subtype_assignment_balanced.csv` | `data/` | Figures 3, 4, 5, 6, 7, 8 |
| `subtype_classifier.json` | `data/` | Figures 3, 5, 7, 8 |
| `Figure_2_consensus_clustering_v7_pubsize.{png,svg}` | `results/` | main text |
| `Supplementary_NMF_validation_v3_pubsize.{png,svg}` | `results/` | supplement |
| `Supplementary_Figure_CDF_analysis_pubsize.{png,svg}` | `results/` | supplement |
| `Supplementary_NMF_validation_metrics_v2.csv` | `results/` | internal |
| `Table1_baseline_characteristics.csv` | `results/` | main text |
| `Table_S2_NCC_Classifier_68genes.xlsx` | `results/` | supplement |
| `survival_blind_analysis_v1_stats.csv` | `results/` | internal |
| `survival_blind_analysis_v1_k2_concordance.csv` | `results/` | internal |
| `survival_blind_analysis_v1_assignment.csv` | `results/` | internal |
| `survival_blind_analysis_v1_crosstab.csv` | `results/` | internal |
| `survival_blind_analysis_v1_consensus_k{2,3,4,5,6}.npy` | `results/` | internal |
| `survival_blind_supplementary_figure_v4.{png,svg}` | `results/figures/` | supplement |
| `classifier_validation.png` | `results/` | internal |

## Run

    cd Figure2_clustering/scripts
    python run_upstream.py

Runtime: 20-60 min. The KMeans consensus clustering runs 100 iterations
per k value, and the NMF validation fits 100 seeds; these dominate the
total time.

Options:

    python run_upstream.py --list              # list steps and exit
    python run_upstream.py --start 3           # skip clustering recomputation
    python run_upstream.py --end 5             # stop after CDF plots
    python run_upstream.py --only nmf_validation
    python run_upstream.py --only survival_blind_analysis_v1
    python run_upstream.py --only survival_blind_supplementary_figure_v4

## After the run

    # to Figure3
    cp data/subtype_assignment_balanced.csv  ../Figure3_transcriptomic_CPTAC/data/processed/
    cp data/subtype_classifier.json          ../Figure3_transcriptomic_CPTAC/data/processed/

    # to Figure4
    cp data/subtype_assignment_balanced.csv  ../Figure4_genomic/data/processed/

    # to Figure5
    cp data/subtype_assignment_balanced.csv  ../Figure5_immune/data/processed/
    cp data/subtype_classifier.json          ../Figure5_immune/data/processed/

    # to Figure6
    cp data/subtype_assignment_balanced.csv  ../Figure6_associations_lipid_immune/data/

    # to Figure7
    cp data/subtype_assignment_balanced.csv  ../Figure7_therapeutic_analysis_pubsize/data/
    cp data/subtype_classifier.json          ../Figure7_therapeutic_analysis_pubsize/data/

    # to Figure8
    cp data/subtype_assignment_balanced.csv  ../Figure8_External_validation_OS_KM_score_combined/data/classifier/
    cp data/subtype_assignment_balanced.csv  ../Figure8_External_validation_OS_KM_score_combined/data/
    cp data/subtype_classifier.json          ../Figure8_External_validation_OS_KM_score_combined/data/classifier/
    cp data/subtype_classifier.json          ../Figure8_External_validation_OS_KM_score_combined/data/

On Windows, use `copy` instead of `cp`.

## Special scripts

- `compute_clustering.py` writes `subtype_assignment_balanced.csv` on first
  run, using the "smaller cluster = Subtype_1" rule. On subsequent runs it
  aligns to the existing file, so it is safe to re-run.
- `survival_blind_analysis_v1.py` is a reviewer-response analysis that
  clusters without any survival information and reports concordance with
  the published subtypes. Its outputs feed
  `survival_blind_supplementary_figure_v4.py`.
- `table1_baseline.py` and `table_s2_ncc_classifier.py` produce table
  outputs (CSV and XLSX) that feed the manuscript's main and supplementary
  tables.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).