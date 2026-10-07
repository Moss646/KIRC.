#!/usr/bin/env python3
# download_gdc_counts.py
# Download GDC STAR counts for TCGA-KIRC primary tumours, merge into a
# gene x sample count matrix, and write the two limma-ready files:
#   KIRC_counts_for_limma.csv
#   sample_info_for_limma.csv
#
# Input:
#   subtype_assignment_balanced.csv  (default: <repo>/data/processed/;
#                                     override with --subtype-csv)
# Output (default: <repo>/data/tcga_kirc/):
#   KIRC_counts_for_limma.csv
#   sample_info_for_limma.csv
#   _gdc_cache/                       per-file tsv cache
#
# Runtime: 10-30 min on first run, seconds when the cache is warm.

import argparse
import io
import json
import os
import sys
import tarfile
import time

import numpy as np
import pandas as pd
import requests

if sys.platform == "win32":
    import io as _io
    sys.stdout = _io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                   errors="replace")


def _find_project_root(start):
    cur = start
    for _ in range(5):
        if os.path.isdir(os.path.join(cur, "data")):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return start


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = _find_project_root(SCRIPT_DIR)

DEFAULT_OUTPUT_DIR = os.path.join(PROJ_ROOT, "data", "tcga_kirc")
DEFAULT_SUBTYPE_CSV = os.path.join(PROJ_ROOT, "data", "processed",
                                   "subtype_assignment_balanced.csv")


def parse_args():
    p = argparse.ArgumentParser(
        description="Download and merge TCGA-KIRC STAR counts for limma.")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                   help="directory for the merged matrix and sample info")
    p.add_argument("--subtype-csv", default=DEFAULT_SUBTYPE_CSV,
                   help="subtype assignment table (Patient, Subtype)")
    p.add_argument("--cache-dir", default=None,
                   help="per-file tsv cache (default: <output-dir>/_gdc_cache)")
    return p.parse_args()


def fetch_file_list():
    print("[1/4] Querying GDC file list...", flush=True)
    r = requests.post("https://api.gdc.cancer.gov/files", json={
        "filters": json.dumps({"op": "and", "content": [
            {"op": "in", "content": {
                "field": "cases.project.project_id", "value": ["TCGA-KIRC"]}},
            {"op": "in", "content": {
                "field": "files.data_type",
                "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {
                "field": "files.analysis.workflow_type",
                "value": ["STAR - Counts"]}},
        ]}),
        "fields": "file_id,cases.samples.submitter_id,cases.samples.sample_type",
        "format": "JSON",
        "size": "2000",
    }, timeout=60)

    uuid_map = {}
    for h in r.json()["data"]["hits"]:
        s = h["cases"][0]["samples"][0]
        if "Primary Tumor" in s.get("sample_type", ""):
            uuid_map[h["file_id"]] = s["submitter_id"][:12]

    print(f"      primary tumours: {len(uuid_map)}", flush=True)
    return uuid_map


def download_tsvs(uuid_map, cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    to_get = [u for u in uuid_map
              if not os.path.exists(os.path.join(cache_dir, f"{u}.tsv"))]
    print(f"[2/4] {len(uuid_map) - len(to_get)} cached, "
          f"{len(to_get)} to download", flush=True)

    if not to_get:
        return

    chunk_size = 50
    total_batches = (len(to_get) - 1) // chunk_size + 1
    for i in range(0, len(to_get), chunk_size):
        chunk = to_get[i:i + chunk_size]
        batch = i // chunk_size + 1
        print(f"      batch {batch}/{total_batches} "
              f"({len(chunk)} files)...", flush=True, end=" ")

        try:
            resp = requests.post("https://api.gdc.cancer.gov/data",
                                 json={"ids": chunk}, timeout=300)
            if resp.status_code != 200:
                print(f"HTTP {resp.status_code}", flush=True)
                continue

            with tarfile.open(fileobj=io.BytesIO(resp.content),
                              mode="r:gz") as tar:
                for member in tar.getmembers():
                    if not (member.isfile() and member.name.endswith(".tsv")):
                        continue
                    fid = member.name.split("/")[0]
                    tar.extract(member, cache_dir)
                    extracted = os.path.join(cache_dir, member.name)
                    target = os.path.join(cache_dir, f"{fid}.tsv")
                    if os.path.exists(extracted) and not os.path.exists(target):
                        os.rename(extracted, target)
                    subdir = os.path.dirname(extracted)
                    if os.path.isdir(subdir) and not os.listdir(subdir):
                        os.rmdir(subdir)
            print("ok", flush=True)
        except Exception as e:
            print(f"failed: {e}", flush=True)

        time.sleep(0.5)


def merge_matrix(uuid_map, cache_dir, subtype_csv):
    print("[3/4] Merging count matrix...", flush=True)
    counts = {}
    for fid, pid in uuid_map.items():
        fpath = os.path.join(cache_dir, f"{fid}.tsv")
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath, sep="\t", comment="#", index_col=0)
        if "unstranded" not in df.columns:
            continue
        ser = df["unstranded"]
        ser = ser[ser.index.astype(str).str.startswith("ENSG")]
        counts[pid] = ser.astype(int)

    matrix = pd.DataFrame(counts)
    matrix = matrix.loc[matrix.sum(axis=1) > 0]
    print(f"      merged: {matrix.shape[0]} genes x {matrix.shape[1]} samples",
          flush=True)

    if not os.path.exists(subtype_csv):
        sys.exit(f"subtype file not found: {subtype_csv}\n"
                 f"Run Figure2 first, or pass --subtype-csv explicitly.")

    sub = pd.read_csv(subtype_csv)
    if not {"Patient", "Subtype"}.issubset(sub.columns):
        sys.exit(f"subtype file must have columns Patient, Subtype: {subtype_csv}")
    sub_map = dict(zip(sub["Patient"], sub["Subtype"]))

    common = sorted(set(matrix.columns) & set(sub_map.keys()))
    group = [sub_map[p] for p in common]
    matrix = matrix[common]

    n1 = sum(1 for g in group if g == "Subtype_1")
    n2 = sum(1 for g in group if g == "Subtype_2")
    print(f"      matched to subtypes: S1={n1}, S2={n2}", flush=True)
    return matrix, pd.DataFrame({"Patient": common, "Subtype": group})


def main():
    args = parse_args()
    cache_dir = args.cache_dir or os.path.join(args.output_dir, "_gdc_cache")
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"  output dir   = {args.output_dir}")
    print(f"  subtype csv  = {args.subtype_csv}")

    uuid_map = fetch_file_list()
    if not uuid_map:
        sys.exit("no primary tumour files found")

    download_tsvs(uuid_map, cache_dir)
    matrix, info = merge_matrix(uuid_map, cache_dir, args.subtype_csv)

    counts_out = os.path.join(args.output_dir, "KIRC_counts_for_limma.csv")
    info_out = os.path.join(args.output_dir, "sample_info_for_limma.csv")
    print("[4/4] Saving...", flush=True)
    matrix.to_csv(counts_out)
    info.to_csv(info_out, index=False)
    print(f"      {counts_out}")
    print(f"      {info_out}", flush=True)


if __name__ == "__main__":
    main()