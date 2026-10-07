"""Path configuration for the Figure 6 reproduction package.

Large inputs (TCGA-KIRC expression matrix) are resolved from the
MEL_DATA_ROOT environment variable, which must point to the directory
holding tcga_kirc/KIRC_expr_log2_tpm.csv. Small bundled inputs live in
../data; outputs are written to ../results.
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

MEL_DATA_ROOT = os.environ.get('MEL_DATA_ROOT', '')

# Large input: TCGA-KIRC expression (log2 TPM+1), genes x samples
EXPR = (os.path.join(MEL_DATA_ROOT, 'tcga_kirc', 'KIRC_expr_log2_tpm.csv')
        if MEL_DATA_ROOT else '')

# Bundled inputs
SUBTYPE_CSV = os.path.join(DATA_DIR, 'subtype_assignment_balanced.csv')
TIDE_CSV = os.path.join(DATA_DIR, 'TIDE_official_results.csv')
XCELL_CSV = os.path.join(DATA_DIR, 'xcell_scores_raw.csv')

# Output of metab_immune_corr_v3.py, input of fig_metab_immune_corr_supp_v4.py
CORR_CSV = os.path.join(RESULTS_DIR, 'metabolic_immune_corr_v2_full.csv')


def require_large_inputs():
    if not (MEL_DATA_ROOT and os.path.isfile(EXPR)):
        raise SystemExit(
            'MEL_DATA_ROOT is not set or tcga_kirc/KIRC_expr_log2_tpm.csv '
            'is missing. Set MEL_DATA_ROOT to the directory containing '
            'tcga_kirc/KIRC_expr_log2_tpm.csv.'
        )