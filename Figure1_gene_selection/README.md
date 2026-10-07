# Figure1_gene_selection

Lipid-metabolism gene selection using MSigDB gene sets, univariate Cox
screening, strict multi-criteria filtering, and balanced layer quotas.
Produces the 68-gene panel used for subtyping in Figure2.

## Inputs

| File | Source | How to get |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | TCGA-KIRC_prep | copy to `data/tcga_kirc/` |
| `KIRC_clinical_cleaned.csv` | TCGA-KIRC_prep | copy to `data/tcga_kirc/` |

MSigDB GMT files and the GO OBO ontology are downloaded automatically
during pipeline steps 1 and 2. Requires internet access.

## Filtering pipeline

1. Download the Hallmark, KEGG, Reactome, and GO:BP collections from
   MSigDB v2026.1.Hs. Keep gene sets whose name or description matches a
   curated lipid keyword list.
2. Deduplicate the C5:GO:BP lipid sets hierarchically using GO ontology
   depth and information content (IC). Remove parents that are covered
   by children and trim by set size.
3. Univariate Cox proportional hazards screening of every lipid gene
   against TCGA-KIRC overall survival; BH FDR.
4. Strict multi-criteria selection:
   p < 0.001, |log2HR| > 0.30, expression SD above median,
   pairwise |r| < 0.85.
5. Balanced 68-gene panel: fixed quota per functional layer
   (Core_Lipid_Enzyme 20, Signaling_Lipoprotein 15,
   Transcriptional_Regulator 15, Lipid_Droplet_Storage 5,
   Ferroptosis_ROS 2, Other 11).

## Outputs

| File | Location | Consumer |
|---|---|---|
| `balanced_genes.txt` | `data/consensus_cluster/` | Figure2, Figure8 |
| `selected_genes_strict.csv` | `data/gene_selection/` | Figure5 |
| `msigdb_lipid_genesets_dedup.csv` | `data/` | Figure3 |
| `msigdb_lipid_gene_pool_dedup.csv` | `data/` | Figure2, Figure8 |
| `msigdb_lipid_metadata.json` | `data/` | internal |
| `go_dedup_stats.json` | `data/` | internal |
| `cox_results_all.csv` | `data/cox_screening/` | internal |
| `cox_results_significant.csv` | `data/cox_screening/` | internal |
| `cox_results_fdr_significant.csv` | `data/cox_screening/` | internal |
| `cox_summary.json` | `data/cox_screening/` | internal |
| `Figure_1_gene_selection_v12_pubsize.{png,svg}` | `results/` | main text |

## Run

    cd Figure1_gene_selection/scripts
    python run_upstream.py              # full pipeline (steps 1-5 + figure)
    python run_upstream.py --start 3    # skip MSigDB / GO downloads
    python run_upstream.py --end 4      # stop before the balanced panel
    python run_upstream.py --no-figure  # skip figure + classifier check

Runtime: 10-20 min, dominated by the univariate Cox loop.

## After the run

    # to Figure2
    cp data/consensus_cluster/balanced_genes.txt      ../Figure2_clustering/data/
    cp data/msigdb_lipid_gene_pool_dedup.csv          ../Figure2_clustering/data/

    # to Figure3
    cp data/msigdb_lipid_genesets_dedup.csv           ../Figure3_transcriptomic_CPTAC/data/gene_sets/

    # to Figure5
    cp data/gene_selection/selected_genes_strict.csv  ../Figure5_immune/data/processed/

    # to Figure8
    cp data/consensus_cluster/balanced_genes.txt      ../Figure8_External_validation_OS_KM_score_combined/data/classifier/
    cp data/msigdb_lipid_gene_pool_dedup.csv          ../Figure8_External_validation_OS_KM_score_combined/data/

On Windows, use `copy` instead of `cp`.

## Citation

If you use this project, please cite the main manuscript. A
machine-readable citation file is available in the repository root
(`CITATION.cff`).