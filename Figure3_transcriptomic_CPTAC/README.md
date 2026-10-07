# Figure3_transcriptomic_CPTAC

Transcriptomic (DEG, GSEA, GSVA) and proteomic (CPTAC) characterization
of the two subtypes. Panel layout: volcano, PCA + PERMANOVA, GSEA,
GSVA, TCGA key gene boxplots, CPTAC protein boxplots, 16-gene
phospholipid remodeling boxplots.

## Inputs

| File | Source | How to get |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/raw/` |
| `KIRC_counts_for_limma.csv` | TCGA-KIRC_prep | copy to `data/raw/` |
| `sample_info_for_limma.csv` | TCGA-KIRC_prep | copy to `data/raw/` |
| `subtype_assignment_balanced.csv` | Figure2_clustering | copy to `data/processed/` |
| `subtype_classifier.json` | Figure2_clustering | copy to `data/processed/` |
| `msigdb_lipid_genesets_dedup.csv` | Figure1_gene_selection | copy to `data/gene_sets/` |
| `rnaseq_tumor.txt` | CPTAC | download to `data/raw/` |
| `proteome_tumor.txt` | CPTAC | download to `data/raw/` |
| `proteome_normal.txt` | CPTAC | download to `data/raw/` |

CPTAC files are available from the CPTAC data portal
(https://proteomic.datacommons.cancer.gov/) or the LinkedOmics portal.

## Pipeline steps

1. `02_differential_exp.R` - edgeR + limma-voom differential expression
   (S1 vs S2). Uses `filterByExpr`, TMM normalization, voom, eBayes.
2. `01_convert_ensembl_to_symbol.R` - map Ensembl gene IDs to HGNC
   symbols using `org.Hs.eg.db`.
3. `make_msigdb_gmt.py` - build the lipid gene set GMT file
   (`msigdb_lipid_dedup_236.gmt`) from the Figure1 deduplicated CSV.
4. `03_fgsea.R` - GSEA on the ranked limma-voom log2FC list.
5. `gsva_R.R` - GSVA pathway scores (Gaussian KCDF).
6. `gsva_limma.R` - limma eBayes on the GSVA scores.
7. `Fig3_main_v18.py` - assemble the 3x3 main figure.
8. `cptac_lipid_remodeling_check_v2.py` - CPTAC protein validation of
   the 16 phospholipid remodeling genes.

## Outputs

| File | Location | Consumer |
|---|---|---|
| `deg_limma_voom.csv` | `data/processed/` | internal |
| `deg_limma_voom_with_symbols.csv` | `data/processed/` | Figures 5, 8 |
| `gsea_full_results.csv` | `data/processed/` | internal |
| `gsva_scores.csv` | `data/processed/` | internal |
| `gsva_subtype_statistics.csv` | `data/processed/` | internal |
| `Figure_3_transcriptomic_CPTAC_v18.{png,svg}` | `results/figures/` | main text |
| `cptac_lipid_remodeling_16genes_check_v2.csv` | `results/tables/` | supplement |

## Run

    cd Figure3_transcriptomic_CPTAC/scripts
    python run_all.py

Runtime: 10-20 min. The `run_all.py` script will set `MEL_DATA_ROOT`
to the project root automatically, so the R scripts resolve paths
relative to the project.

## After the run

    # to Figure5
    cp data/processed/deg_limma_voom_with_symbols.csv  ../Figure5_immune/data/processed/

    # to Figure8
    cp data/processed/deg_limma_voom_with_symbols.csv  ../Figure8_External_validation_OS_KM_score_combined/data/deg/
    cp data/processed/deg_limma_voom_with_symbols.csv  ../Figure8_External_validation_OS_KM_score_combined/data/

On Windows, use `copy` instead of `cp`.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).