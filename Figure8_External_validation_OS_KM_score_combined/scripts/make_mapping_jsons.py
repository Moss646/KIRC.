#!/usr/bin/env python3
# make_mapping_jsons.py
# Generate the two gene -> Ensembl ID mapping JSON files required by the
# Figure8 external validation scripts.
#
# Input (auto-discovered):
#   Figure1's balanced_genes.txt   68-gene NCC classifier
# Output:
#   data/gene_ensembl_mapping.json
#   data/lipid16_ensembl_mapping_v1.json
#
# Requires internet access to rest.ensembl.org for symbols not covered by
# the bundled static table. Use --offline to skip network calls entirely.

import argparse
import json
import os
import sys
import time
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(PROJ_ROOT)


# ----------------------------- gene panels -----------------------------

KEY14 = [
    "CPT1A", "ACOX1", "CPT2", "ACADSB", "ACADM", "ACAA2", "CD36",
    "SLC27A2", "FASN", "SCD", "PLIN2", "PPARG", "PPARGC1A", "ITPKA",
]

LIPID16 = [
    "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "MBOAT1", "MBOAT2",
    "LCLAT1", "MBOAT7", "PLA2G4A", "PLA2G6", "PNPLA8", "PLA2G4C",
    "PLA2G2A", "PLA2G4F", "PLA2G1B", "PLA2G2D",
]

NON_NCC9 = [
    "CPT1A", "ACOX1", "CPT2", "CD36", "SLC27A2",
    "FASN", "SCD", "PLIN2", "PPARG",
]

# Static Ensembl IDs for everything except the 68 NCC genes coming from
# Figure1. These are stable across Ensembl releases.
STATIC = {
    # key14
    "CPT1A":    "ENSG00000110090",
    "ACOX1":    "ENSG00000161533",
    "CPT2":     "ENSG00000157184",
    "ACADSB":   "ENSG00000196177",
    "ACADM":    "ENSG00000117054",
    "ACAA2":    "ENSG00000167315",
    "CD36":     "ENSG00000135218",
    "SLC27A2":  "ENSG00000140284",
    "FASN":     "ENSG00000169710",
    "SCD":      "ENSG00000099194",
    "PLIN2":    "ENSG00000147872",
    "PPARG":    "ENSG00000132170",
    "PPARGC1A": "ENSG00000109819",
    "ITPKA":    "ENSG00000137825",
    # lipid16 (only the ones not already listed above)
    "LPCAT1":   "ENSG00000153395",
    "LPCAT2":   "ENSG00000042813",
    "LPCAT3":   "ENSG00000111684",
    "LPCAT4":   "ENSG00000175893",
    "MBOAT1":   "ENSG00000172197",
    "MBOAT2":   "ENSG00000143797",
    "LCLAT1":   "ENSG00000172954",
    "MBOAT7":   "ENSG00000125505",
    "PLA2G4A":  "ENSG00000116711",
    "PLA2G6":   "ENSG00000184381",
    "PNPLA8":   "ENSG00000135241",
    "PLA2G4C":  "ENSG00000105499",
    "PLA2G2A":  "ENSG00000188257",
    "PLA2G4F":  "ENSG00000168907",
    "PLA2G1B":  "ENSG00000170890",
    "PLA2G2D":  "ENSG00000117207",
}


# ----------------------------- helpers -----------------------------

def find_balanced_genes():
    """Locate balanced_genes.txt in this repo or in a sibling Figure1 folder."""
    candidates = [
        os.path.join(PROJ_ROOT, "data", "balanced_genes.txt"),
        os.path.join(PROJ_ROOT, "data", "classifier", "balanced_genes.txt"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c

    for d in sorted(os.listdir(REPO_ROOT)):
        if not d.startswith("Figure1"):
            continue
        for dirpath, _, files in os.walk(os.path.join(REPO_ROOT, d)):
            if "balanced_genes.txt" in files:
                return os.path.join(dirpath, "balanced_genes.txt")
    return None


def ensembl_lookup(symbol, retries=3, timeout=10):
    url = (f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/"
           f"{symbol}?content-type=application/json")
    for _ in range(retries):
        try:
            req = urllib.request.Request(url,
                                         headers={"User-Agent": "mapper"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = json.loads(r.read().decode())
                eid = data.get("id")
                if eid:
                    return eid
        except Exception:
            time.sleep(1.0)
    return None


def build_mapping(genes, allow_network):
    out = {}
    to_query = []

    for g in sorted(set(genes)):
        if g in STATIC:
            out[g] = STATIC[g]
        else:
            to_query.append(g)

    if to_query and allow_network:
        print(f"  querying Ensembl for {len(to_query)} symbols "
              f"({len(out)} resolved from static table)")
        for i, g in enumerate(to_query, 1):
            eid = ensembl_lookup(g)
            if eid:
                out[g] = eid
                print(f"    [{i}/{len(to_query)}] {g:<12s} -> {eid}")
            else:
                print(f"    [{i}/{len(to_query)}] {g:<12s} -> NOT FOUND")
            time.sleep(0.05)
    elif to_query:
        print(f"  [offline] {len(to_query)} symbols not covered by the "
              f"static table: {to_query}")

    return out


# ----------------------------- main -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Generate gene -> Ensembl mapping JSONs for Figure8.")
    ap.add_argument("--out-dir", default=os.path.join(PROJ_ROOT, "data"),
                    help="directory for the two output JSON files")
    ap.add_argument("--offline", action="store_true",
                    help="do not call Ensembl; use bundled static table only")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    balanced = find_balanced_genes()
    ncc68 = []
    if balanced:
        with open(balanced) as f:
            ncc68 = [line.strip() for line in f if line.strip()]
        print(f"[1/3] NCC genes from {balanced}: {len(ncc68)}")
    else:
        print("[1/3] balanced_genes.txt not found; "
              "run Figure1 first for full coverage")

    all_genes = set(ncc68) | set(KEY14) | set(LIPID16) | set(NON_NCC9)
    mapping = build_mapping(all_genes, allow_network=not args.offline)

    out1 = os.path.join(args.out_dir, "gene_ensembl_mapping.json")
    with open(out1, "w", encoding="utf-8") as f:
        json.dump(mapping, f, indent=2, sort_keys=True)
    print(f"[2/3] wrote {out1} ({len(mapping)} genes)")

    lipid = {g: mapping[g] for g in LIPID16 if g in mapping}
    out2 = os.path.join(args.out_dir, "lipid16_ensembl_mapping_v1.json")
    with open(out2, "w", encoding="utf-8") as f:
        json.dump(lipid, f, indent=2, sort_keys=True)
    print(f"[3/3] wrote {out2} ({len(lipid)}/16 genes)")

    if len(lipid) < 16:
        missing = [g for g in LIPID16 if g not in lipid]
        print(f"  WARNING: missing lipid16 entries: {missing}")


if __name__ == "__main__":
    main()