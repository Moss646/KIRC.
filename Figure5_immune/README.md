# Figure5_immune

Immune microenvironment analysis: CIBERSORT, xCell, ssGSEA, GSVA,
camera, purity-stratified marker genes, HIF/hypoxia sensitivity
analyses, and scRNA-seq cell-type expression.

## Inputs

### Project inputs (copy from other projects)

| File | Source | How to get |
|---|---|---|
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/processed/` |
| `subtype_classifier.json` | Figure2_clustering | copy to `data/processed/` |
| `selected_genes_strict.csv` | Figure1_gene_selection | copy to `data/processed/` |
| `deg_limma_voom_with_symbols.csv` | Figure3_transcriptomic_CPTAC | copy to `data/processed/` |

### Upstream inputs (from TCGA-KIRC_prep)

| File | Source | How to get |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/tcga_kirc/` |
| `KIRC_expr_tpm.csv` | TCGA-KIRC_prep | copy to `data/tcga_kirc/` |
| `KIRC_counts_for_limma.csv` | TCGA-KIRC_prep | copy to `data/raw/` |
| `sample_info_for_limma.csv` | TCGA-KIRC_prep | copy to `data/raw/` |

### External downloads

| File | Source | How to get |
|---|---|---|
| `h.all.v2024.1.Hs.symbols.gmt` | MSigDB | download to `data/gene_sets/` |
| `c2.all.v2024.1.Hs.symbols.gmt` | MSigDB | download to `data/gene_sets/` |
| `c5.all.v2023.2.Hs.symbols.gmt` | MSigDB | download to `data/gene_sets/` |
| `GSE159115_*.h5` | GEO | download to `data/scrnaseq/` |
| `GSE159115_ccRCC_anno.csv.gz` | GEO | download to `data/scrnaseq/` |

The C5 GO:BP gmt file must be from release 2023.2 - the "response to
hypoxia" GO terms were retired from the ontology in late 2023 and are
absent from v2024.1.

## Outputs

### R / Python intermediate

| File | Location | Consumer |
|---|---|---|
| `cibersort_results_official_full.csv` | `data/processed/` | internal (Panel A of `fig4_immune_merged_v17.py`) |
| `cibersort_results_official_filtered.csv` | `data/processed/` | internal |
| `cibersort_results.csv` | `data/processed/` | internal |
| `xcell_scores_raw.csv` | `data/xcell/` | Figure6 |
| `xcell_scores_zscore.csv` | `data/xcell/` | internal |
| `xcell_subtype_diff.csv` | `data/xcell/` | internal |
| `xcell_spill.rds` | `data/xcell/` | internal |
| `ssgsea_immune_scores_R.csv` | `data/processed/` | internal |
| `ssgsea_immune_stats_R.csv` | `data/processed/` | internal |
| `gsva_hif_sets_scores_v1.csv` | `data/processed/` | internal |
| `ssgsea_hif_sets_scores_v1.csv` | `data/processed/` | internal |
| `immune_purity_strata_v3.csv` | `data/processed/` | internal |
| `Supplementary_Table_NCC_68_genes.xlsx` | `data/processed/` | Figure5 internal |
| `sample_subtype_mapping_v5.csv` | `results/tables/` | internal |

### Figures and tables

| File | Location | Consumer |
|---|---|---|
| `fig4_immune_merged_v17.{png,svg}` | `results/figures/` | supplement |
| `hif_gene_sets_three_methods_v1.{png,svg}` | `results/figures/` | supplement |
| `hif_pathway_activity_v2.{png,svg}` | `results/figures/` | supplement |
| `Supplementary_immune_features_composition_independence_v4.{png,svg}` | `results/figures/` | supplement |
| `Supplementary_xCell_analysis_v17_pubsize.{png,svg}` | `results/figures/` | supplement |
| `Fig_supp_scrnaseq_v5_pubsize.{png,svg}` | `results/figures/` | supplement |
| `lipid_remodeling_16genes_S1S2_v5.{png,svg}` | `results/figures/` | supplement |
| `pla2_inflammatory_S1S2_v1.{png,svg}` | `results/figures/` | supplement |
| `hif_gene_sets_stats_v1.csv` | `results/tables/` | internal |
| `hif_gene_sets_gsva_stats_v1.csv` | `results/tables/` | internal |
| `hif_gene_sets_camera_stats_v1.csv` | `results/tables/` | internal |
| `hif_gene_sets_zscore_stats_v1.csv` | `results/tables/` | internal |
| `hif_gene_sets_per_gene_check_v1.csv` | `results/tables/` | internal |
| `immune_features_composition_stratified_mwu_v2.csv` | `results/tables/` | internal |
| `immune_genes_stratified_limma_v3.csv` | `results/tables/` | internal |
| `immune_genes_voom_log2cpm_all_v3.csv` | `data/processed/` | internal |
| `lipid_remodeling_16genes_stats_v5.csv` | `results/tables/` | internal |
| `pla2_inflammatory_stats_v1.csv` | `results/tables/` | internal |
| `fig4_immune_merged_C_D_limma_v16.csv` | `results/tables/` | internal |
| `panel68_cell_source_binning_v1.csv` | `results/tables/` | internal |
| `v5_checkpoint_subtype.csv` | `results/tables/` | internal |
| `v5_mhc_subtype.csv` | `results/tables/` | internal |

## Run

    cd Figure5_immune/scripts
    python run_all.py

Runtime: 3-8 h. The three scRNA-seq scripts that read the GSE159115 h5
files dominate this time. If their outputs already exist, you can delete
those three scripts from the `script_order` list to rerun the rest in
minutes.

The first script (`extract_spill.R`) only needs to be run once to
extract xCell's internal spillover matrix to an RDS file. It will be
skipped automatically on subsequent runs.

## After the run

    # to Figure6
    cp data/xcell/xcell_scores_raw.csv  ../Figure6_associations_lipid_immune/data/

On Windows, use `copy` instead of `cp`.

## Special scripts

- `run_cibersort_official.R` runs CIBERSORT with `perm=1000` and no
  random seed. Re-running it changes the P-value cutoff and therefore
  the CIBERSORT sample subset. If you want reproducible downstream
  results, add `set.seed(42)` before the `cibersort(...)` call.
- `01_pseudobulk_ncc_v5.py`, `02_celltype_subtype_expression_v5.py`, and
  `panel68_cell_source_binning_v1.py` each read every h5 file once.
  They are the three slowest scripts in this project.
- `make_purity_strata_v3.py` and `make_stratified_mwu_v2.py` include a
  hash check that writes a `_recomputed.csv` file if the new result
  differs from the existing CSV. Inspect the diff before replacing the
  canonical file.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).