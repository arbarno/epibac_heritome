#!/usr/bin/env python3

"""
> calc_mean_methylation.py <

Multi-functional script calculates, on a per-gene basis, two values:
1. methylation density (dense vs. sparse methylation)
2. mean methylation %

Requires a genome, a GFF3 file with gene information, and an annotated
bed file (to extract # of methylated CpGs in each gene).
This is a split script from calc_bias_density_medians.py
"""

import argparse
import csv
import statistics
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="Calculate methylation density and mean methylation percentage.")
parser.add_argument('genome', metavar='FASTA file', type=argparse.FileType('r'), help='FASTA file of genome.')
parser.add_argument('gff', metavar='GFF3 file', type=argparse.FileType('r'), help='GFF3 file of genome.')
parser.add_argument('bed', metavar='Annotated BED file', type=argparse.FileType('r'), help='BED file with methylated positions.')
args = parser.parse_args()

# Parse genome, GFF3 file, and BED file
genome_sequence = parse_fasta.get_all_sequences(args.genome, 'fasta')
scaffold_gff3 = parse_gff3.parse_gff3(args.gff, 'gene')

bed_data = {}
tsv_reader = csv.reader(args.bed, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    gene_id = row[6]
    meth_pct = float(row[4])
    if gene_id not in bed_data:
        bed_data[gene_id] = []
    bed_data[gene_id].append(meth_pct)

print('gene', 'meth pos', 'meth density', 'mean meth %', sep='\t')
for scaf in scaffold_gff3:
    for gene in scaffold_gff3[scaf]:
        if gene not in bed_data:
            print(gene, 0, 0, 0, sep='\t')
        else:
            meth_density = round(len(bed_data[gene]) / 2 * 100, 2)
            mean_meth_pct = round(statistics.mean(bed_data[gene]), 2)
            print(gene, len(bed_data[gene]), meth_density, mean_meth_pct, sep='\t')
