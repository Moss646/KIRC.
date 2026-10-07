"""Path configuration. All paths resolved relative to this file.

Layout assumed:
    <repo>/scripts/config.py
    <repo>/data/
    <repo>/results/
"""

import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

CLASSIFIER_DIR = DATA_DIR / "classifier"
CPTAC_DIR = DATA_DIR / "cptac"
GEO_DIR = DATA_DIR / "geo"
ICGC_DIR = DATA_DIR / "icgc"
DEG_DIR = DATA_DIR / "deg"


def _pick(*candidates):
    for p in candidates:
        if Path(p).exists():
            return str(p)
    return str(candidates[0])


SUBTYPE_CLASSIFIER_JSON = _pick(
    CLASSIFIER_DIR / "subtype_classifier.json",
    DATA_DIR / "subtype_classifier.json",
)
BALANCED_GENES_TXT = _pick(
    CLASSIFIER_DIR / "balanced_genes.txt",
    DATA_DIR / "balanced_genes.txt",
)
SUBTYPE_ASSIGNMENT_CSV = _pick(
    CLASSIFIER_DIR / "subtype_assignment_balanced.csv",
    DATA_DIR / "subtype_assignment_balanced.csv",
)
TCGA_CLINICAL_CSV = _pick(
    CLASSIFIER_DIR / "KIRC_clinical_cleaned.csv",
    DATA_DIR / "KIRC_clinical_cleaned.csv",
)
TCGA_DEG_CSV = _pick(
    DEG_DIR / "deg_limma_voom_with_symbols.csv",
    DATA_DIR / "deg_limma_voom_with_symbols.csv",
)

TCGA_EXPR = _pick(
    DATA_DIR / "tcga_kirc" / "KIRC_expr_log2_tpm.csv",
    DATA_DIR / "KIRC_expr_log2_tpm.csv",
)

CPTAC_RNA = _pick(
    DATA_DIR / "rnaseq_tumor.txt",
    CPTAC_DIR / "rnaseq_tumor.txt",
)
CPTAC_PROTEIN = _pick(
    CPTAC_DIR / "proteome_tumor.txt",
    DATA_DIR / "proteome_tumor.txt",
)
CPTAC_CLINICAL = _pick(
    CPTAC_DIR / "clinical.txt",
    DATA_DIR / "clinical.txt",
)

GEO_EMTAB_EXPR = _pick(
    GEO_DIR / "EMTAB1980_expr.txt",
    DATA_DIR / "EMTAB1980_expr.txt",
)
GEO_NM2GENE = _pick(
    GEO_DIR / "emtrab_nm_to_gene.json",
    DATA_DIR / "emtrab_nm_to_gene.json",
)
GEO_SATO_CLINICAL = _pick(
    GEO_DIR / "Sato2013_TableS1.xlsx",
    DATA_DIR / "Sato2013_TableS1.xlsx",
)

ICGC_EXP_SEQ = _pick(
    ICGC_DIR / "exp_seq.RECA-EU.tsv.gz",
    DATA_DIR / "exp_seq.RECA-EU.tsv.gz",
)
ICGC_SPECIMEN = _pick(
    ICGC_DIR / "specimen.RECA-EU.tsv.gz",
    DATA_DIR / "specimen.RECA-EU.tsv.gz",
)
ICGC_DONOR = _pick(
    ICGC_DIR / "donor.RECA-EU.tsv.gz",
    DATA_DIR / "donor.RECA-EU.tsv.gz",
)

GENE_MAPPING_JSON = _pick(
    DATA_DIR / "gene_ensembl_mapping.json",
    DATA_DIR / "gene_mapping.json",
)
LIPID16_ENSEMBL_JSON = _pick(
    DATA_DIR / "lipid16_ensembl_mapping_v1.json",
    DATA_DIR / "lipid16_ensembl_mapping.json",
)

RESULTS_DIR = str(REPO_ROOT / "results")
TABLES_DIR = str(REPO_ROOT / "results" / "tables")
FIGURES_DIR = str(REPO_ROOT / "results" / "figures")
PROVENANCE_DIR = str(REPO_ROOT / "results" / "tables" / "_provenance")
WORK = TABLES_DIR

for _d in (TABLES_DIR, FIGURES_DIR, PROVENANCE_DIR):
    os.makedirs(_d, exist_ok=True)


def _data_version():
    try:
        h = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        if h.returncode == 0 and h.stdout.strip():
            return h.stdout.strip()
    except Exception:
        pass
    return os.environ.get("NCC_DATA_VERSION", "unknown")


DATA_VERSION = _data_version()


def emit_provenance(script_name, rows):
    payload = {
        "script": script_name,
        "data_version": DATA_VERSION,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "rows": rows,
    }
    out = os.path.join(PROVENANCE_DIR, Path(script_name).stem + ".json")
    try:
        with open(out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        print(f"[PROVENANCE] wrote {len(rows)} rows -> {out}")
    except Exception as e:
        print(f"[PROVENANCE] WARNING: failed to write {out}: {e}")


def get_gene_ensembl_mapping():
    with open(GENE_MAPPING_JSON) as f:
        return json.load(f)


def print_paths():
    print("NCC external validation - paths")
    print(f"  REPO_ROOT      = {REPO_ROOT}")
    print(f"  DATA_DIR       = {DATA_DIR}")
    print(f"  CLASSIFIER     = {SUBTYPE_CLASSIFIER_JSON}")
    print(f"  TCGA_EXPR      = {TCGA_EXPR}")
    print(f"  CPTAC_RNA      = {CPTAC_RNA}")
    print(f"  CPTAC_PROTEIN  = {CPTAC_PROTEIN}")
    print(f"  GEO_EMTAB_EXPR = {GEO_EMTAB_EXPR}")
    print(f"  ICGC_EXP_SEQ   = {ICGC_EXP_SEQ}")
    print(f"  GENE_MAPPING   = {GENE_MAPPING_JSON}")
    print(f"  TABLES_DIR     = {TABLES_DIR}")