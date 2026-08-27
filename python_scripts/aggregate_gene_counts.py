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


def process_sample(bed_path: str) -> dict:
    """Returns {gene: {'meth': int, 'total': int, 'n_cpgs': int}}
    Matches the logic of the original calc_mean_methylation.py: gene_id is
    taken directly from column 7 (row[6]). For genic rows this is the real
    gene ID; for intergenic rows this is literally the string "intergenic",
    which will not match any real gene ID and is excluded naturally by the
    downstream gene_universe filter"""
    gene_data = defaultdict(lambda: {"meth": 0, "total": 0, "n_cpgs": 0})

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

            gene_data[gene_id]["meth"] += n_mod
            gene_data[gene_id]["total"] += n_valid_cov
            gene_data[gene_id]["n_cpgs"] += 1

    return gene_data


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-CpG counts to per-gene methylated/total counts."
    )
    parser.add_argument("bed_files", nargs="+",
                         help="Annotated bed files, one per sample "
                              "(e.g. *.10x.clair3.filt.conservedCpGs.annotated.bed)")
    parser.add_argument("--out-tsv", default="gene_level_counts.tsv")
    args = parser.parse_args()

    with open(args.out_tsv, "w") as out:
        out.write("sample\tgene\tn_cpgs\tmethylated_count\ttotal_count\tweighted_average\n")

        for bed_path in args.bed_files:
            sample_name = os.path.basename(bed_path)
            # strip common suffix to get a clean sample name
            sample_name = re.sub(
                r"\.10x\.clair3\.filt\.conservedCpGs\.annotated\.bed$", "", sample_name
            )

            print(f"Processing {sample_name} ({bed_path})...")
            gene_data = process_sample(bed_path)

            for gene, counts in gene_data.items():
                if counts["total"] > 0:
                    weighted_avg = round(100 * counts["meth"] / counts["total"], 4)
                else:
                    weighted_avg = 0.0

                out.write(f"{sample_name}\t{gene}\t{counts['n_cpgs']}\t"
                          f"{counts['meth']}\t{counts['total']}\t{weighted_avg}\n")

    print(f"\nDone. Output written to: {args.out_tsv}")


if __name__ == "__main__":
    main()
