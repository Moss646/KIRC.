#!/usr/bin/env python3
"""Download lipid-related MSigDB gene sets and build a gene pool."""

import csv
import io
import json
import sys
import time
from collections import Counter
from pathlib import Path

import requests

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
OUTPUT_DIR = REPO_ROOT / "data"

MSIGDB_BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs"

COLLECTIONS = {
    "Hallmark": "h.all.v2026.1.Hs.symbols.gmt",
    "C2:KEGG": "c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt",
    "C2:REACTOME": "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
    "C5:GO:BP": "c5.go.bp.v2026.1.Hs.symbols.gmt",
}

LIPID_KEYWORDS = [
    "LIPID", "FATTY_ACID", "CHOLESTEROL", "TRIGLYCERIDE",
    "PHOSPHOLIPID", "SPHINGOLIPID", "LIPOPROTEIN", "PEROXISOME",
    "ADIPOGENESIS", "STEROID_BIOSYNTHESIS", "ACYL", "LIPASE",
    "CARNITINE", "KETONE_BODY", "BILE_ACID", "EICOSANOID",
    "GLYCEROLIPID", "GLYCEROPHOSPHOLIPID", "ARACHIDONIC",
    "LEUKOTRIENE", "CERAMIDE", "DIACYLGLYCEROL", "CHYLOMICRON",
    "HDL_PARTICLE", "LDL", "VLDL", "LIPID_DROPLET",
    "LIPID_STORAGE", "LIPID_OXIDATION", "FATTY_ACID_BETA",
    "FATTY_ACID_OMEGA", "FATTY_ACID_ALPHA",
    "LIPID_TRANSPORT", "LIPID_LOCALIZATION",
    "PHOSPHATIDYL", "LYSOPHOSPHOLIPID", "PROSTAGLANDIN",
    "TERPENOID", "STEROID_HORMONE",
    "BILE_ACID_AND_BILE_SALT",
    "LIPID_CATABOLIC", "LIPID_BIOSYNTHETIC",
    "FERROPTOSIS",
    "PPAR_SIGNALING",
]


def download_gmt(url, filename):
    full_url = f"{MSIGDB_BASE}/{url}"
    local_path = OUTPUT_DIR / filename

    print(f"downloading {filename} ... ", end="", flush=True)
    resp = requests.get(full_url, timeout=120)
    if resp.status_code != 200:
        print(f"failed (HTTP {resp.status_code})")
        return None

    with open(local_path, "w", encoding="utf-8") as f:
        f.write(resp.text)
    print(f"ok ({len(resp.text):,} bytes)")
    return local_path


def parse_gmt(filepath):
    genesets = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            description = parts[1] if len(parts) > 1 else ""
            genes = [g.strip() for g in parts[2:] if g.strip()]
            genesets.append({
                "name": name,
                "description": description,
                "genes": genes,
                "num_genes": len(genes),
            })
    return genesets


def match_lipid_keywords(name, description):
    text = (name + " " + description).upper()
    return any(kw.upper() in text for kw in LIPID_KEYWORDS)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("MSigDB lipid gene-set download")
    downloaded = {}
    for label, gmt_file in COLLECTIONS.items():
        local = download_gmt(gmt_file, gmt_file)
        if local:
            downloaded[label] = local
        else:
            print(f"skip {label}")

    if not downloaded:
        raise RuntimeError("all downloads failed")

    all_genesets = []
    collection_stats = {}

    for label, filepath in downloaded.items():
        genesets = parse_gmt(filepath)
        lipid_matched = [
            gs for gs in genesets
            if match_lipid_keywords(gs["name"], gs["description"])
        ]
        for gs in lipid_matched:
            gs["source"] = label

        all_genesets.extend(lipid_matched)
        collection_stats[label] = {
            "total": len(genesets),
            "lipid": len(lipid_matched),
        }
        print(f"{label}: {len(genesets)} sets, matched {len(lipid_matched)}")

    all_genes = set()
    gene_source_count = Counter()
    gene_source_list = {}

    for gs in all_genesets:
        for gene in gs["genes"]:
            all_genes.add(gene)
            gene_source_count[gene] += 1
            gene_source_list.setdefault(gene, []).append(gs["name"])

    total_before_dedup = sum(gs["num_genes"] for gs in all_genesets)
    multi_source = sum(1 for c in gene_source_count.values() if c >= 2)
    multi_3plus = sum(1 for c in gene_source_count.values() if c >= 3)

    print(f"genes before dedup: {total_before_dedup}")
    print(f"unique genes: {len(all_genes)}")
    print(f"genes in >=2 sets: {multi_source}")
    print(f"genes in >=3 sets: {multi_3plus}")

    summary_file = OUTPUT_DIR / "msigdb_lipid_genesets_summary.csv"
    with open(summary_file, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["Source", "Geneset_Name", "Num_Genes", "Description", "Genes"])
        for gs in sorted(all_genesets, key=lambda x: (x["source"], x["name"])):
            writer.writerow([
                gs["source"],
                gs["name"],
                gs["num_genes"],
                gs["description"][:300],
                "|".join(gs["genes"]),
            ])

    for label in ["Hallmark", "C2:KEGG", "C2:REACTOME", "C5:GO:BP"]:
        gsets = [gs for gs in all_genesets if gs["source"] == label]
        if not gsets:
            continue
        list_file = OUTPUT_DIR / (
            f"msigdb_{label.replace(':', '_').replace('.', '_')}_lipidsets.txt"
        )
        with open(list_file, "w", encoding="utf-8") as f:
            f.write(f"# {label} lipid-related gene sets\n")
            f.write(f"# total: {len(gsets)}\n")
            f.write("# MSigDB v2026.1.Hs\n\n")
            for gs in sorted(gsets, key=lambda x: x["name"]):
                f.write(f"# {gs['name']} (n={gs['num_genes']})\n")
                f.write("\t".join(gs["genes"]) + "\n\n")

    gene_file = OUTPUT_DIR / "msigdb_lipid_gene_pool.csv"
    with open(gene_file, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene_Symbol", "Num_GeneSets", "Source_GeneSets"])
        for gene in sorted(all_genes):
            writer.writerow([
                gene,
                gene_source_count[gene],
                "; ".join(gene_source_list[gene][:10]),
            ])

    gmt_file = OUTPUT_DIR / "msigdb_lipid_merged.gmt"
    with open(gmt_file, "w", encoding="utf-8") as f:
        for gs in sorted(all_genesets, key=lambda x: (x["source"], x["name"])):
            clean_name = f"{gs['source']}__{gs['name']}"
            f.write(f"{clean_name}\t{gs['description'][:100]}\t")
            f.write("\t".join(gs["genes"]) + "\n")

    metadata = {
        "msigdb_version": "2026.1.Hs",
        "download_date": time.strftime("%Y-%m-%d"),
        "total_genesets": len(all_genesets),
        "total_unique_genes": len(all_genes),
        "collection_summary": {
            label: {
                "total_genesets": collection_stats[label]["total"],
                "lipid_genesets": collection_stats[label]["lipid"],
                "unique_genes": len(set(
                    g for gs in all_genesets if gs["source"] == label
                    for g in gs["genes"]
                )),
            }
            for label in collection_stats
        },
        "genesets": [
            {
                "source": gs["source"],
                "name": gs["name"],
                "num_genes": gs["num_genes"],
                "description": gs["description"][:200],
            }
            for gs in sorted(all_genesets, key=lambda x: (x["source"], x["name"]))
        ],
        "top_frequent_genes": [
            {"gene": g, "count": c}
            for g, c in gene_source_count.most_common(50)
        ],
        "gene_pool": sorted(all_genes),
    }

    json_file = OUTPUT_DIR / "msigdb_lipid_metadata.json"
    with open(json_file, "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print(f"saved {summary_file}")
    print(f"saved {gene_file}")
    print(f"saved {gmt_file}")
    print(f"saved {json_file}")


if __name__ == "__main__":
    main()