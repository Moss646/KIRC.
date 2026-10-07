#!/usr/bin/env python3
"""Hierarchical deduplication of C5:GO:BP lipid gene sets."""

import csv
import gzip
import io
import json
import math
import sys
from collections import defaultdict, deque
from pathlib import Path

import requests

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
OUTPUT_DIR = REPO_ROOT / "data"
GENESETS_CSV = OUTPUT_DIR / "msigdb_lipid_genesets_summary.csv"
GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

MIN_IC = 1.2
MAX_IC = 5.5
MIN_GENES = 15
MAX_GENES = 500
MAX_ANCESTOR_COVERAGE = 0.55


def parse_go_obo(text):
    terms = {}
    current = None

    for line in text.split("\n"):
        line = line.rstrip("\r")

        if line == "[Term]":
            if current and current.get("id"):
                terms[current["id"]] = current
            current = {"id": "", "name": "", "namespace": "", "parents": []}
        elif line == "[Typedef]":
            if current and current.get("id"):
                terms[current["id"]] = current
            current = None
        elif current is not None:
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:]
            elif line.startswith("is_a: "):
                parent = line[6:].split("!")[0].strip()
                current["parents"].append(parent)

    if current and current.get("id"):
        terms[current["id"]] = current

    bp = {go_id: t for go_id, t in terms.items() if t.get("namespace") == "biological_process"}
    print(f"parsed GO OBO: {len(terms)} terms, {len(bp)} BP terms")
    return bp


def build_name_index(bp_terms):
    idx = {}
    for go_id, t in bp_terms.items():
        key = t["name"].upper().replace(" ", "_").replace("-", "_")
        idx[key] = go_id
    return idx


def load_gobp_genesets():
    genesets = {}
    with open(GENESETS_CSV, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["Source"] != "C5:GO:BP":
                continue
            genes = [g.strip() for g in row["Genes"].split("|") if g.strip()]
            genesets[row["Geneset_Name"]] = set(genes)
    return genesets


def match_gobp_to_go(genesets, bp_terms):
    name_to_id = build_name_index(bp_terms)
    matched = {}
    unmatched = []

    for gs_name in genesets:
        simple = gs_name.replace("GOBP_", "")
        if simple in name_to_id:
            matched[gs_name] = name_to_id[simple]
        else:
            unmatched.append(gs_name)

    print(f"mapped: {len(matched)}/{len(genesets)}")
    if unmatched:
        print(f"unmapped: {len(unmatched)}")
        for u in unmatched[:10]:
            print(f"  {u}")
    return matched


def get_ancestors(go_id, bp_terms, max_depth=50):
    ancestors = set()
    queue = deque([go_id])
    depth = 0

    while queue and depth < max_depth:
        for _ in range(len(queue)):
            cur = queue.popleft()
            if cur in bp_terms:
                for parent in bp_terms[cur].get("parents", []):
                    if parent not in ancestors:
                        ancestors.add(parent)
                        queue.append(parent)
        depth += 1

    return ancestors


def dedup(go_matched, genesets, bp_terms):
    total_bp_genes = len(set().union(*genesets.values()))
    info = {}

    for gs_name, go_id in go_matched.items():
        genes = genesets[gs_name]
        n = len(genes)
        ic = -math.log10(max(n, 1) / total_bp_genes) if n > 0 else 0
        info[gs_name] = {
            "go_id": go_id,
            "go_name": bp_terms.get(go_id, {}).get("name", gs_name),
            "n_genes": n,
            "ic": round(ic, 2),
            "genes": genes,
            "ancestors": get_ancestors(go_id, bp_terms),
            "keep": True,
            "reason": "",
        }

    go_to_gs = {v["go_id"]: k for k, v in info.items()}
    chain_map = defaultdict(set)

    for gs_name, d in info.items():
        root = gs_name
        for anc_id in sorted(
            d["ancestors"],
            key=lambda a: len(get_ancestors(a, bp_terms)),
            reverse=True,
        ):
            if anc_id in go_to_gs:
                root = go_to_gs[anc_id]
        chain_map[root].add(gs_name)

    removed_count = 0
    for root, members in chain_map.items():
        if len(members) < 2:
            continue

        chain = sorted(members, key=lambda n: info[n]["ic"])
        for i in range(len(chain)):
            if not info[chain[i]]["keep"]:
                continue

            parent_genes = info[chain[i]]["genes"]
            if not parent_genes:
                continue

            for j in range(i + 1, len(chain)):
                if not info[chain[j]]["keep"]:
                    continue

                child_genes = info[chain[j]]["genes"]
                coverage = len(parent_genes & child_genes) / len(parent_genes)

                if coverage > MAX_ANCESTOR_COVERAGE:
                    info[chain[i]]["keep"] = False
                    info[chain[i]]["reason"] = (
                        f"covered by {info[chain[j]]['go_name'][:50]} "
                        f"({coverage:.0%})"
                    )
                    removed_count += 1
                    break

    for gs_name, d in info.items():
        if not d["keep"]:
            continue
        if d["ic"] < MIN_IC:
            d["keep"] = False
            d["reason"] = f"IC={d['ic']:.2f} < {MIN_IC}"
            removed_count += 1
        elif d["ic"] > MAX_IC:
            d["keep"] = False
            d["reason"] = f"IC={d['ic']:.2f} > {MAX_IC}"
            removed_count += 1
        elif d["n_genes"] < MIN_GENES:
            d["keep"] = False
            d["reason"] = f"n={d['n_genes']} < {MIN_GENES}"
            removed_count += 1
        elif d["n_genes"] > MAX_GENES:
            d["keep"] = False
            d["reason"] = f"n={d['n_genes']} > {MAX_GENES}"
            removed_count += 1

    kept = {k: v for k, v in info.items() if v["keep"]}
    removed = {k: v for k, v in info.items() if not v["keep"]}
    print(f"kept {len(kept)}, removed {len(removed)}")
    return kept, removed, info


def save_results(kept, removed, all_info):
    with open(OUTPUT_DIR / "go_bp_lipid_dedup_kept.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["GOBP_Name", "GO_ID", "GO_Name", "Num_Genes", "IC", "Genes"])
        for k in sorted(kept, key=lambda x: kept[x]["ic"]):
            v = kept[k]
            w.writerow([k, v["go_id"], v["go_name"], v["n_genes"], v["ic"], "|".join(sorted(v["genes"]))])

    with open(OUTPUT_DIR / "go_bp_lipid_dedup_removed.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["GOBP_Name", "GO_ID", "GO_Name", "Num_Genes", "IC", "Reason"])
        for k in sorted(removed, key=lambda x: removed[x]["ic"]):
            v = removed[k]
            w.writerow([k, v["go_id"], v["go_name"], v["n_genes"], v["ic"], v["reason"]])

    with open(OUTPUT_DIR / "go_bp_lipid_dedup_full_report.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["GOBP_Name", "GO_ID", "GO_Name", "Num_Genes", "IC", "Status", "Reason"])
        for k, v in sorted(all_info.items(), key=lambda x: x[1]["ic"]):
            w.writerow([k, v["go_id"], v["go_name"], v["n_genes"], v["ic"], "KEPT" if v["keep"] else "REMOVED", v["reason"]])

    merged_rows = []
    with open(GENESETS_CSV, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["Source"] != "C5:GO:BP":
                merged_rows.append(row)

    for k, v in kept.items():
        merged_rows.append({
            "Source": "C5:GO:BP",
            "Geneset_Name": k,
            "Num_Genes": str(v["n_genes"]),
            "Description": v["go_name"],
            "Genes": "|".join(sorted(v["genes"])),
        })

    merged_file = OUTPUT_DIR / "msigdb_lipid_genesets_dedup.csv"
    with open(merged_file, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Source", "Geneset_Name", "Num_Genes", "Description", "Genes"])
        for row in merged_rows:
            w.writerow([
                row["Source"],
                row["Geneset_Name"],
                row["Num_Genes"],
                row["Description"][:300],
                row["Genes"],
            ])

    all_genes = set()
    for row in merged_rows:
        for g in row["Genes"].split("|"):
            g = g.strip()
            if g:
                all_genes.add(g)

    pool_file = OUTPUT_DIR / "msigdb_lipid_gene_pool_dedup.csv"
    with open(pool_file, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Gene_Symbol"])
        for g in sorted(all_genes):
            w.writerow([g])

    stats = {
        "original_gobp_count": len(all_info),
        "kept": len(kept),
        "removed": len(removed),
        "final_total_genesets": len(merged_rows),
        "final_gene_pool": len(all_genes),
    }
    with open(OUTPUT_DIR / "go_dedup_stats.json", "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)

    print(f"saved kept: {len(kept)}")
    print(f"saved removed: {len(removed)}")
    print(f"saved gene pool: {len(all_genes)}")


def main():
    print("GO:BP hierarchical deduplication")

    try:
        resp = requests.get(GO_OBO_URL + ".gz", timeout=120)
        if resp.status_code == 200:
            text = gzip.decompress(resp.content).decode("utf-8")
        else:
            raise RuntimeError("gzip download failed")
    except Exception:
        resp = requests.get(GO_OBO_URL, timeout=120)
        text = resp.text

    bp_terms = parse_go_obo(text)
    genesets = load_gobp_genesets()
    print(f"loaded {len(genesets)} GO:BP lipid gene sets")

    go_matched = match_gobp_to_go(genesets, bp_terms)
    if not go_matched:
        raise RuntimeError("no GO terms matched")

    kept, removed, all_info = dedup(go_matched, genesets, bp_terms)
    save_results(kept, removed, all_info)

    print(f"done: {len(all_info)} -> {len(kept)}")


if __name__ == "__main__":
    main()