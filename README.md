# KIRC Lipid-Metabolism Subtypes

A multi-cohort bioinformatics study identifying and validating two
lipid-metabolism subtypes (S1 and S2) of clear cell renal cell carcinoma
(ccRCC / TCGA-KIRC), with validation in E-MTAB-1980, ICGC RECA-EU, and
CPTAC cohorts.

## Repository structure

| Folder | Purpose | Entry point |
|---|---|---|
| `TCGA-KIRC_prep/` | Raw TCGA-KIRC cleaning + GDC raw counts | `tcga_kirc_prep.py`, `download_gdc_counts.py` |
| `Figure1_gene_selection/` | Gene selection (MSigDB + Cox + strict filter) | `scripts/run_upstream.py` |
| `Figure2_clustering/` | Consensus clustering + NCC classifier | `scripts/run_upstream.py` |
| `Figure3_transcriptomic_CPTAC/` | DEG / GSEA / GSVA / CPTAC | `scripts/run_all.py` |
| `Figure4_genomic/` | Mutation and CNA landscape | `scripts/run_all.py` |
| `Figure5_immune/` | CIBERSORT / xCell / ssGSEA / scRNA-seq | `scripts/run_all.py` |
| `Figure6_associations_lipid_immune/` | Lipid-immune correlation | `scripts/run_all.py` |
| `Figure7_therapeutic_analysis_pubsize/` | oncoPredict + TIDE | `scripts/run_all.py` |
| `Figure8_External_validation_OS_KM_score_combined/` | Multi-cohort external validation | `scripts/run_all.py` |

## Dependency graph
TCGA-KIRC_prep
│
├──→ Figure1 ──→ balanced_genes.txt ──→ Figure2
│ │
│ ├──→ selected_genes_strict.csv ──→ Figure5
│ │
│ └──→ msigdb_lipid_genesets_dedup.csv ──→ Figure3
│
├──→ Figure2 ──→ subtype_assignment_balanced.csv,
│ subtype_classifier.json
│ │
│ ├──→ Figure3 ──→ deg_limma_voom_with_symbols.csv ──→ Figure5, Figure8
│ ├──→ Figure4
│ ├──→ Figure5 ──→ xcell_scores_raw.csv ──→ Figure6
│ ├──→ Figure6
│ ├──→ Figure7 ──→ TIDE_official_results.csv ──→ Figure6
│ └──→ Figure8 ──→ nonncc9_fao_score_per_sample_v1.csv ──→ Figure4

## Pipeline execution

Run projects in this order. After each project finishes, copy its
outputs to the downstream projects listed below.

### 1. TCGA-KIRC_prep

```bash
cd TCGA-KIRC_prep
python tcga_kirc_prep.py
python download_gdc_counts.py --subtype-csv ../Figure2_clustering/data/subtype_assignment_balanced.csv
copy data/tcga_kirc/KIRC_expr_log2_tpm.csv to Figure1, Figure2,
Figure3, Figure5, Figure6, Figure7, Figure8

copy data/tcga_kirc/KIRC_expr_tpm.csv to Figure5, Figure8

copy data/tcga_kirc/KIRC_clinical_cleaned.csv to Figure2, Figure3,
Figure8

copy data/tcga_kirc/KIRC_counts_for_limma.csv and
data/tcga_kirc/sample_info_for_limma.csv to Figure3 (data/raw/) and
Figure5 (data/raw/)

download_gdc_counts.py must be run after Figure2 has produced
subtype_assignment_balanced.csv.

2. Figure1_gene_selection
bash
cd Figure1_gene_selection/scripts
python run_upstream.py
copy data/consensus_cluster/balanced_genes.txt to
Figure2 (data/) and Figure8 (data/classifier/)

copy data/gene_selection/selected_genes_strict.csv to Figure5
(data/processed/)

copy data/msigdb_lipid_genesets_dedup.csv to Figure3
(data/gene_sets/)

copy data/msigdb_lipid_gene_pool_dedup.csv to Figure2 (data/) and
Figure8 (data/)

3. Figure2_clustering
bash
cd Figure2_clustering/scripts
python run_upstream.py
copy data/subtype_assignment_balanced.csv to Figure3
(data/processed/), Figure4 (data/processed/), Figure5
(data/processed/), Figure6 (data/), Figure7 (data/), Figure8
(data/classifier/ and data/)

copy data/subtype_classifier.json to Figure3 (data/processed/),
Figure5 (data/processed/), Figure7 (data/), Figure8
(data/classifier/ and data/)

4. Figure3_transcriptomic_CPTAC
bash
cd Figure3_transcriptomic_CPTAC/scripts
python run_all.py
copy data/processed/deg_limma_voom_with_symbols.csv to Figure5
(data/processed/) and Figure8 (data/deg/ and data/)

5. Figure5_immune
bash
cd Figure5_immune/scripts
python run_all.py
copy data/xcell/xcell_scores_raw.csv to Figure6 (data/)

6. Figure7_therapeutic_analysis_pubsize
bash
cd Figure7_therapeutic_analysis_pubsize/scripts
python run_all.py
copy data/processed/TIDE_official_results.csv to Figure6 (data/)

7. Figure8_External_validation_OS_KM_score_combined
bash
cd Figure8_External_validation_OS_KM_score_combined/scripts
python run_all.py
The first step (make_mapping_jsons.py) generates two gene → Ensembl
mapping JSON files inside data/. It requires balanced_genes.txt from
Figure1.

copy results/tables/nonncc9_fao_score_per_sample_v1.csv to
Figure4 (data/raw/)

8. Figure6_associations_lipid_immune
bash
cd Figure6_associations_lipid_immune/scripts
python run_all.py
9. Figure4_genomic
bash
cd Figure4_genomic/scripts
python run_all.py
If a downstream script reports a missing file, check the "Inputs"
section of its README to see which upstream project produces it and
where it should be copied.

Environment
Python 3.10+ and R 4.4+ are both required.

bash
python -m pip install -r requirements.txt
Rscript install_packages.R
Data
Large raw files (TCGA downloads, GDSC2 RDS, scRNA-seq h5, external
cohorts) are not included. Each project's README lists its required
inputs and their sources.

Citation
See CITATION.cff.

License
MIT. See LICENSE.