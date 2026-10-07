# TCGA-KIRC_prep

Upstream cleaning of raw TCGA-KIRC downloads. Produces the standard
expression matrix, cleaned clinical table, and limma-ready raw counts
shared by all downstream figures.

## Inputs

### External downloads (from GDC)

| File | Source | How to get |
|---|---|---|
| `TCGA-KIRC_gene_annotation.csv` | GDC | https://portal.gdc.cancer.gov/ |
| `mrna/TCGA-KIRC_tpm_mrna.csv` | GDC | https://portal.gdc.cancer.gov/ |
| `TCGA-KIRC_clinical_survival.csv` | GDC | https://portal.gdc.cancer.gov/ |
| `TCGA-KIRC_clinical_info.csv` | GDC | https://portal.gdc.cancer.gov/ |

Place the four files under `data/raw/`. Note that `TCGA-KIRC_tpm_mrna.csv`
must be inside a `mrna/` subdirectory.

### Project input

| File | Source | How to get |
|---|---|---|
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/processed/` |

`subtype_assignment_balanced.csv` is only required by
`download_gdc_counts.py`, and only after Figure2 has produced it.

## Cleaning steps applied by `tcga_kirc_prep.py`

1. Keep protein-coding genes only.
2. Keep primary tumour samples only (`sample_code == '01'`).
3. Keep the earliest vial per patient (`01A` before `01B`).
4. Low-expression filter: mean TPM >= 1.0 AND expressed in >= 20 % of samples.

Step 4 changes the gene count and is easy to omit from Methods descriptions
in a manuscript; note it when reporting the final matrix dimensions.

## Outputs

| File | Location | Consumer |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | `data/tcga_kirc/` | Figures 1, 2, 3, 5, 6, 7, 8 |
| `KIRC_expr_tpm.csv` | `data/tcga_kirc/` | Figures 5, 8 |
| `KIRC_clinical_cleaned.csv` | `data/tcga_kirc/` | Figures 2, 3, 8 |
| `KIRC_data_summary.json` | `data/tcga_kirc/` | internal |
| `KIRC_counts_for_limma.csv` | `data/tcga_kirc/` | Figures 3, 5 |
| `sample_info_for_limma.csv` | `data/tcga_kirc/` | Figures 3, 5 |

## Run

    cd TCGA-KIRC_prep

    # 1. Clean expression + clinical (no network access required)
    python tcga_kirc_prep.py

    # 2. Download GDC STAR raw counts (requires Figure2 output and internet)
    python download_gdc_counts.py --subtype-csv ../Figure2_clustering/data/subtype_assignment_balanced.csv

If Figure2 has not yet been run, copy `subtype_assignment_balanced.csv` into
`data/processed/` first, then run `download_gdc_counts.py` without the
`--subtype-csv` flag.

`download_gdc_counts.py` caches per-file TSVs under
`data/tcga_kirc/_gdc_cache/`. Delete the cache to force a fresh download.

## After the run

Copy the outputs to their downstream consumers:

    # expression matrices
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure1_gene_selection/data/tcga_kirc/
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure2_clustering/data/
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure3_transcriptomic_CPTAC/data/raw/
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure5_immune/data/tcga_kirc/
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure6_associations_lipid_immune/data/tcga_kirc/
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure7_therapeutic_analysis_pubsize/data/
    cp data/tcga_kirc/KIRC_expr_log2_tpm.csv      ../Figure8_External_validation_OS_KM_score_combined/data/tcga_kirc/

    cp data/tcga_kirc/KIRC_expr_tpm.csv           ../Figure5_immune/data/tcga_kirc/
    cp data/tcga_kirc/KIRC_expr_tpm.csv           ../Figure8_External_validation_OS_KM_score_combined/data/tcga_kirc/

    # clinical table
    cp data/tcga_kirc/KIRC_clinical_cleaned.csv   ../Figure2_clustering/data/
    cp data/tcga_kirc/KIRC_clinical_cleaned.csv   ../Figure3_transcriptomic_CPTAC/data/processed/
    cp data/tcga_kirc/KIRC_clinical_cleaned.csv   ../Figure8_External_validation_OS_KM_score_combined/data/classifier/

    # raw counts for limma
    cp data/tcga_kirc/KIRC_counts_for_limma.csv   ../Figure3_transcriptomic_CPTAC/data/raw/
    cp data/tcga_kirc/KIRC_counts_for_limma.csv   ../Figure5_immune/data/raw/
    cp data/tcga_kirc/sample_info_for_limma.csv   ../Figure3_transcriptomic_CPTAC/data/raw/
    cp data/tcga_kirc/sample_info_for_limma.csv   ../Figure5_immune/data/raw/

On Windows, use `copy` instead of `cp`.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).