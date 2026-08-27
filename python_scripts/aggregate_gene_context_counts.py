#!/usr/bin/env python3
"""
aggregate_gene_counts.py

Purpose
-------
Aggregate per-CpG methylation counts (Nmod, Nvalid_cov) to per-gene totals,
for each sample, to produce the methylated/total count pairs needed for a
proper binomial/quasi-binomial GLM (rather than a pre-averaged percentage).

Input
-----
One or more *.conservedCpGs.annotated.bed files (one per sample), format:
  col1: scaffold
  col2: start
  col3: end
  col4: Nvalid_cov (total coverage)
  col5: percent methylated
  col6: Nmod (methylated count)
  col7: genomic context
  col8-11: distance/position annotation
  col12: gene annotation, gene ID after '|' (e.g. "downstream_watson|g1")

Output
------
One TSV with columns: sample, gene, methylated_count, total_count, n_cpgs
(n_cpgs = number of CpG positions contributing to this gene's totals -
equivalent to the old 'meth_pos' column, kept for consistency with the
>=3 CpG gene-inclusion filter used elsewhere in the pipeline)
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict


def extract_gene_id(annotation_field: str) -> str:
    """Extract gene ID from a field like 'downstream_watson|g1' -> 'g1'."""
    if "|" in annotation_field:
        return annotation_field.split("|")[-1].strip()
    return annotation_field.strip()


def classify_context(context_field: str, sub_annotation_field: str = None) -> str:
    """Classify a row's genomic context as intergenic, exon, intron, or other.
    - Intergenic rows: column 7 (row[6]) is literally "intergenic"
    - Genic rows: column 7 is a real gene ID; sub-context (Exon_N_.../Intron_N_...)
      is in column 11 (row[10])"""
    if context_field.strip() == "intergenic":
        return "intergenic"
    if sub_annotation_field is None:
        return "other"
    if sub_annotation_field.startswith("Exon"):
        return "exon"
    if sub_annotation_field.startswith("Intron"):
        return "intron"
    return "other"


def process_sample(bed_path: str) -> tuple:
    """Returns (gene_data, context_data):
    gene_data: {gene: {'meth': int, 'total': int, 'n_cpgs': int}}
    context_data: {context: {'meth': int, 'total': int, 'n_cpgs': int}}
    Both built in a single pass over the file."""
    gene_data = defaultdict(lambda: {"meth": 0, "total": 0, "n_cpgs": 0})
    context_data = defaultdict(lambda: {"meth": 0, "total": 0, "n_cpgs": 0})

    with open(bed_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 7:
                continue
            try:
                n_valid_cov = int(row[3])
                n_mod = int(row[5])
                gene_id = row[6].strip()
            except (ValueError, IndexError):
                continue

            if gene_id.startswith("tRNA"):
                continue

            # gene-level aggregation (unchanged - only real gene IDs, not "intergenic")
            if gene_id != "intergenic":
                gene_data[gene_id]["meth"] += n_mod
                gene_data[gene_id]["total"] += n_valid_cov
                gene_data[gene_id]["n_cpgs"] += 1

            # context-level aggregation (new)
            sub_annotation = row[10].strip() if len(row) > 10 else None
            context = classify_context(gene_id, sub_annotation)
            context_data[context]["meth"] += n_mod
            context_data[context]["total"] += n_valid_cov
            context_data[context]["n_cpgs"] += 1

    # Add a "total" category: sum of intergenic + exon + intron only,
    # explicitly excluding "other" (ambiguous/no_info positions)
    real_contexts = ["intergenic", "exon", "intron"]
    total_meth = sum(context_data[c]["meth"] for c in real_contexts if c in context_data)
    total_cov = sum(context_data[c]["total"] for c in real_contexts if c in context_data)
    total_n_cpgs = sum(context_data[c]["n_cpgs"] for c in real_contexts if c in context_data)
    context_data["total"] = {"meth": total_meth, "total": total_cov, "n_cpgs": total_n_cpgs}

    return gene_data, context_data


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-CpG counts to per-gene and per-context methylated/total counts."
    )
    parser.add_argument("bed_files", nargs="+",
                         help="Annotated bed files, one per sample "
                              "(e.g. *.10x.clair3.filt.conservedCpGs.annotated.bed)")
    parser.add_argument("--out-gene-tsv", default="gene_level_counts.tsv")
    parser.add_argument("--out-context-tsv", default="context_level_counts.tsv")
    args = parser.parse_args()

    gene_out = open(args.out_gene_tsv, "w")
    context_out = open(args.out_context_tsv, "w")

    gene_out.write("sample\tgene\tn_cpgs\tmethylated_count\ttotal_count\tweighted_average\n")
    context_out.write("sample\tcontext\tn_cpgs\tmethylated_count\ttotal_count\tweighted_average\n")

    for bed_path in args.bed_files:
        sample_name = os.path.basename(bed_path)
        sample_name = re.sub(
            r"\.10x\.clair3\.filt\.conservedCpGs\.annotated\.bed$", "", sample_name
        )

        print(f"Processing {sample_name} ({bed_path})...")
        gene_data, context_data = process_sample(bed_path)

        for gene, counts in gene_data.items():
            weighted_avg = round(100 * counts["meth"] / counts["total"], 4) if counts["total"] > 0 else 0.0
            gene_out.write(f"{sample_name}\t{gene}\t{counts['n_cpgs']}\t"
                            f"{counts['meth']}\t{counts['total']}\t{weighted_avg}\n")

        for context, counts in context_data.items():
            weighted_avg = round(100 * counts["meth"] / counts["total"], 4) if counts["total"] > 0 else 0.0
            context_out.write(f"{sample_name}\t{context}\t{counts['n_cpgs']}\t"
                               f"{counts['meth']}\t{counts['total']}\t{weighted_avg}\n")

    gene_out.close()
    context_out.close()

    print(f"\nDone.")
    print(f"Gene-level output written to: {args.out_gene_tsv}")
    print(f"Context-level output written to: {args.out_context_tsv}")


if __name__ == "__main__":
    main()
