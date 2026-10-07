"""Build fao9_immune_input_tcga_v1.csv from the Figure8 FAO9 score table.

Reads the per-sample FAO9 score table that was copied from Figure8 into
data/raw/, keeps the TCGA rows only, and writes the two-column table that
driver_mut_fao_v1.py and fig_driver_3p_del_FAO9_v1.py expect:

    data/processed/fao9_immune_input_tcga_v1.csv    (sample, FAO9)

Run this once before the Figure4 driver-mutation scripts.
"""

import os
import sys

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)

SRC = os.path.join(PROJ_ROOT, "data", "raw",
                   "nonncc9_fao_score_per_sample_v1.csv")
OUT_DIR = os.path.join(PROJ_ROOT, "data", "processed")
OUT = os.path.join(OUT_DIR, "fao9_immune_input_tcga_v1.csv")


def main():
    if not os.path.exists(SRC):
        sys.exit(f"missing input: {SRC}\n"
                 f"Copy nonncc9_fao_score_per_sample_v1.csv from\n"
                 f"  Figure8_External_validation_OS_KM_score_combined/"
                 f"results/tables/\n"
                 f"into Figure4_genomic/data/raw/ before running this script.")

    df = pd.read_csv(SRC)
    required = {"cohort", "sample", "score"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"input missing columns: {sorted(missing)}")

    tcga = df[df["cohort"] == "TCGA"][["sample", "score"]].copy()
    tcga = tcga.rename(columns={"score": "FAO9"})
    tcga = tcga.sort_values("sample").reset_index(drop=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    tcga.to_csv(OUT, index=False)

    print(f"wrote {OUT} ({len(tcga)} TCGA samples)")


if __name__ == "__main__":
    main()