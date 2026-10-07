\# Data



Raw and intermediate data are \*\*not tracked\*\*. Download the files below

and place them under the matching sub-directories before running the

pipeline.



\## `data/tcga\_kirc/`



| File | Source |

|---|---|

| `KIRC\_expr\_log2\_tpm.csv` | \[TCGA GDC](https://portal.gdc.cancer.gov/) project TCGA-KIRC, log2(TPM+1) |

| `KIRC\_expr\_tpm.csv` | Same source, raw TPM |



\## `data/processed/`



| File | Source |

|---|---|

| `subtype\_assignment\_balanced.csv` | NCC subtype labels (`Patient`, `Subtype`) |

| `deg\_limma\_voom\_with\_symbols.csv` | From `Figure3\_transcriptomic\_CPTAC` |

| `Supplementary\_Table\_NCC\_68\_genes.xlsx` | From `Figure1\_gene\_selection` |

| `cibersort\_results\_official\_full.csv` | CIBERSORT output (spillover-corrected) |

| `ssgsea\_immune\_scores\_R.csv`, `ssgsea\_immune\_stats\_R.csv` | ssGSEA scores and stats from R |

| `immune\_genes\_stratified\_limma\_v3.csv` | Stratified limma output |



\## `data/xcell/`



| File | Source |

|---|---|

| `xcell\_scores\_raw.csv` | R xCell v1.1.0 raw scores |

| `xcell\_scores\_zscore.csv` | Row-wise z-scores |

| `xcell\_subtype\_diff.csv` | S1 vs S2 differential table |



\## `data/scrnaseq/`



GSE159115 ccRCC scRNA-seq:



\- `GSE159115\_ccRCC\_anno.csv.gz` — cell annotations

\- `GSM\*\_SI\_\*.h5` — one 10x HDF5 file per sample (7 samples)



Download from \[GEO GSE159115](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159115).



\## `data/gene\_sets/`



MSigDB gene sets (symbols):



\- `h.all.v2024.1.Hs.symbols.gmt`

\- `c2.all.v2024.1.Hs.symbols.gmt`

\- `c5.all.v2023.2.Hs.symbols.gmt`



Download from \[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/).



\## Sibling-repo dependencies



The following two files are produced by other repositories in this

project. Copy them into `data/processed/` before running:



```bash

copy ..\\Figure3\_transcriptomic\_CPTAC\\data\\processed\\deg\_limma\_voom\_with\_symbols.csv data\\processed\\

copy ..\\Figure1\_gene\_selection\\results\\Supplementary\_Table\_NCC\_68\_genes.xlsx data\\processed\\

```



Alternatively, set the environment variables `LIMMA\_CSV` and

`PANEL68\_XLSX` to point at those files in place:



```bat

set LIMMA\_CSV=D:\\github0918\\Figure3\_transcriptomic\_CPTAC\\data\\processed\\deg\_limma\_voom\_with\_symbols.csv

set PANEL68\_XLSX=D:\\github0918\\Figure1\_gene\_selection\\results\\Supplementary\_Table\_NCC\_68\_genes.xlsx

```



The `run\_all.py` runner already sets both of these when spawning child

processes.

