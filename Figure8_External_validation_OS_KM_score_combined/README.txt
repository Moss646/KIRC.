# Figure8_External_validation_OS_KM_score_combined

Multi-cohort external validation of the 68-gene NCC classifier and the
14-gene lipid-metabolism panel across E-MTAB-1980, ICGC RECA-EU, and
CPTAC. Generates prognostic Cox models (univariate, stage-stratified,
bootstrap C-index), KM curves, cross-cohort classifier concordance, and
the 9-gene non-NCC FAO score used by Figure4.

## Inputs

### Project inputs (copy from other projects)

| File | Source | How to get |
|---|---|---|
| `subtype_classifier.json` | Figure2_clustering | copy to `data/classifier/` or `data/` |
| `balanced_genes.txt` | Figure1_gene_selection | copy to `data/classifier/` or `data/` |
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/classifier/` or `data/` |
| `KIRC_clinical_cleaned.csv` | TCGA-KIRC_prep | copy to `data/classifier/` or `data/` |
| `deg_limma_voom_with_symbols.csv` | Figure3_transcriptomic_CPTAC | copy to `data/deg/` or `data/` |
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/tcga_kirc/` or `data/` |
| `msigdb_lipid_gene_pool_dedup.csv` | Figure1_gene_selection | copy to `data/` |

### External downloads

| File | Source | How to get |
|---|---|---|
| `rnaseq_tumor.txt` | CPTAC | download to `data/cptac/` or `data/` |
| `proteome_tumor.txt` | CPTAC | download to `data/cptac/` or `data/` |
| `clinical.txt` | CPTAC | download to `data/cptac/` or `data/` |
| `EMTAB1980_expr.txt` | GEO | download to `data/geo/` or `data/` |
| `emtrab_nm_to_gene.json` | GEO | download to `data/geo/` or `data/` |
| `Sato2013_TableS1.xlsx` | GEO | download to `data/geo/` or `data/` |
| `exp_seq.RECA-EU.tsv.gz` | ICGC | download to `data/icgc/` or `data/` |
| `specimen.RECA-EU.tsv.gz` | ICGC | download to `data/icgc/` or `data/` |
| `donor.RECA-EU.tsv.gz` | ICGC | download to `data/icgc/` or `data/` |
| `TCGA-KIRC.merged.maf` | GDC | download to `data/` |
| `TCGA-KIRC.gistic.tsv` | GDAC Firehose | download to `data/` |

The two gene -> Ensembl mapping JSON files
(`gene_ensembl_mapping.json` and `lipid16_ensembl_mapping_v1.json`) are
generated automatically by `make_mapping_jsons.py` during the PREP step.

## Pipeline stages (`run_all.py`)

**PREP** (1 step)

1. `make_mapping_jsons.py` - build the two mapping JSON files from
   Figure1's `balanced_genes.txt` plus a bundled static table. Requires
   internet access to rest.ensembl.org for the 68 NCC genes.

**CORE** (7 steps)

2. `simple_cox_all_cohorts.py` - simple Cox HR / log-rank P for the
   three cohorts with OS data.
3. `prognostic_cox.py` - stage-stratified Cox and bootstrap C-index,
   one figure per cohort.
4. `external_validation_KM_pubsize.py` - publication-size KM curves.
5. `classifier_concordance.py` - 68-gene Pearson concordance vs TCGA.
6. `key14_biological.py` - 14-gene direction concordance and Pearson
   correlation with the TCGA effect sizes.
7. `panel_G_CPTAC_protein.py` - CPTAC protein validation of the 14-gene
   panel.
8. `diagnostics.py` - PH assumption check and survival-time filter audit.

**EXTRA** (6 steps)

9. `key14_cross_cohort_v2.py` - 14-gene cross-cohort v2 (uniform
   mean-difference reference).
10. `lipid16_cross_cohort_v2.py` - 16-gene lipid remodeling cross-cohort
    v2.
11. `lipid31_cross_cohort_v1.py` - 30-gene combined cross-cohort v1.
12. `nonncc9_fao_cross_cohort_v1.py` - 9-gene non-NCC FAO score across
    four cohorts. **This produces the FAO9 table consumed by Figure4.**
13. `external_validation_KM_score_combined_v3.py` - combined KM + FAO
    score figure.
14. `driver_mut_fao_v1.py` - driver mutations vs FAO score (this script
    is a duplicate of Figure4's `driver_mut_fao_v1.py`; keep it for
    self-containedness, or remove it from EXTRA if you prefer to run
    it only from Figure4).

**SUMMARY** (2 steps)

15. `verify_all_values.py` - cross-check every value in the summary CSV
    against recomputed values.
16. `generate_summary_csv.py` - assemble
    `NCC_multi_cohort_validation_summary.csv` from provenance JSON
    sidecars.

**DOCX** (3 steps)

17. `results/tables/build_table_docx.py` - main Table 1.
18. `results/tables/build_split_tables.py` - Tables 1, S1, S2.
19. `results/tables/build_Table3.py` - Table 3.

## Outputs

| File | Location | Consumer |
|---|---|---|
| `gene_ensembl_mapping.json` | `data/` | internal |
| `lipid16_ensembl_mapping_v1.json` | `data/` | internal |
| `nonncc9_fao_score_per_sample_v1.csv` | `results/tables/` | Figure4 |
| `nonncc9_fao_score_stats_v1.csv` | `results/tables/` | internal |
| `nonncc9_genes_comparison_v1.csv` | `results/tables/` | internal |
| `nonncc9_cross_cohort_summary_v1.csv` | `results/tables/` | internal |
| `classifier68_cross_cohort_summary_v2.csv` | `results/tables/` | internal |
| `classifier68_genes_comparison_v2.csv` | `results/tables/` | internal |
| `key14_genes_comparison.csv` | `results/tables/` | internal |
| `key14_genes_comparison_v2.csv` | `results/tables/` | internal |
| `key14_cross_cohort_summary_v2.csv` | `results/tables/` | internal |
| `lipid16_genes_comparison_v2.csv` | `results/tables/` | internal |
| `lipid16_cross_cohort_summary_v2.csv` | `results/tables/` | internal |
| `lipid31_genes_comparison_v1.csv` | `results/tables/` | internal |
| `lipid31_cross_cohort_summary_v1.csv` | `results/tables/` | internal |
| `driver_mut_fao_v1.csv` | `results/tables/` | internal |
| `driver_gistic_del_v1.csv` | `results/tables/` | internal |
| `NCC_multi_cohort_validation_summary.csv` | `results/tables/` | main table |
| `_provenance/*.json` | `results/tables/_provenance/` | internal |
| `NCC_validation_Table1.docx` | `results/tables/` | main text |
| `Table1_prognostic.docx` | `results/tables/` | main text |
| `TableS1_classifier.docx` | `results/tables/` | supplement |
| `TableS2_biological.docx` | `results/tables/` | supplement |
| `Table3_reproducibility.docx` | `results/tables/` | main text |
| `Figure_7_prognostic_cox_stratified.{png,svg}` | `results/figures/` | main text |
| `EMTAB_prognostic_cox.{png,svg}` | `results/figures/` | main text |
| `ICGC_prognostic_cox.{png,svg}` | `results/figures/` | main text |
| `External_validation_OS_KM_pubsize.{png,svg}` | `results/figures/` | main text |
| `External_validation_OS_KM_score_combined_v3.{png,svg}` | `results/figures/` | main text |
| `nonncc9_fao_score_4cohorts_v1.{png,svg}` | `results/figures/` | main text |
| `Figure_G_CPTAC_protein_v1.{png,svg}` | `results/figures/` | supplement |

## Run

    cd Figure8_External_validation_OS_KM_score_combined/scripts
    python run_all.py

Options:

    python run_all.py --skip-prep        # skip make_mapping_jsons.py (JSONs already exist)
    python run_all.py --skip-docx        # analysis only, skip Word table generation
    python run_all.py --docx-only        # regenerate Word tables only
    python run_all.py --only prognostic_cox.py
    python run_all.py --only make_mapping_jsons.py

Runtime: 30-60 min. The `prognostic_cox.py` step runs 2000 bootstrap
C-index iterations per cohort, which dominates.

## Requirements

- `Rscript` must be on `PATH`:
  `nonncc9_fao_cross_cohort_v1.py` calls
  `nonncc9_limma_external_v1.R` as a subprocess.
- Internet access for the PREP step (Ensembl REST lookups).
- All Python dependencies listed in the repository `requirements.txt`.

## After the run

    # to Figure4
    cp results/tables/nonncc9_fao_score_per_sample_v1.csv  ../Figure4_genomic/data/raw/

On Windows, use `copy` instead of `cp`.

## Special scripts

- `_recompute_table2_pvals.py` is a standalone diagnostic that
  recomputes every Table 2 P value and writes
  `_recompute_table2_pvals.txt` in the scripts directory. It is not part
  of `run_all.py`; run it manually when you want to cross-check the
  table.
- `diagnostics.py` performs a proportional-hazards check (Schoenfeld
  residuals) and audits the survival-time filter across cohorts. It is
  read-only and prints to stdout only.
- `config.py` includes a `DATA_VERSION` field that is filled from the
  current Git commit hash. If you clone the repository without `.git`,
  set the `NCC_DATA_VERSION` environment variable to whatever identifier
  you prefer.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).