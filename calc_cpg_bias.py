#!/usr/bin/env python3

"""
> calc_cpg_bias.py <

This script calculates, on a per-gene basis, three values:
1. the CpG bias (CpG O/E) ratios for individual genes in the genome

Requires a genome (to calculate expected CpG fraction)and a GFF3 file with gene
information (to calculate observed CpG fraction).
This is a split script from calc_bias_density_medians.py
"""

import argparse
import collections
import itertools
import parse_fasta
import parse_gff3

def reverse_complement(seq):
    seq = seq.replace('U', 'T')
    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)
    return seq[::-1].translate(translation_table)

def get_gene_sequence(sequence, coords):
    startpos, endpos = coords
    revcomp_seq = False
    if startpos > endpos:
        startpos, endpos = endpos, startpos
        revcomp_seq = True
    sliced_seq = sequence[startpos:endpos]
    return reverse_complement(sliced_seq) if revcomp_seq else sliced_seq

parser = argparse.ArgumentParser(description="Calculate CpG bias for genes.")
parser.add_argument('genome', metavar='FASTA file', type=argparse.FileType('r'), help='FASTA file of genome.')
parser.add_argument('gff', metavar='GFF3 file', type=argparse.FileType('r'), help='GFF3 file of genome.')
args = parser.parse_args()

# Parse genome and GFF3 file
genome_sequence = parse_fasta.get_all_sequences(args.genome, 'fasta')
scaffold_gff3 = parse_gff3.parse_gff3(args.gff, 'gene')

# Perform per-gene CpG bias calculations
valid_dinuc = [''.join(x) for x in itertools.product('ACGT', repeat=2)]
cx_dinuc = ['C' + x for x in 'ACGT']
gx_dinuc = ['G' + x for x in 'ACGT']

print('gene', 'observed CpG', 'expected CpG', 'CpG bias', sep='\t')
for scaf in scaffold_gff3:
    for gene in scaffold_gff3[scaf]:
        gene_sequence = get_gene_sequence(genome_sequence[scaf], scaffold_gff3[scaf][gene].coords)
        gene_dinuc = collections.Counter([''.join(x) for x in zip(gene_sequence[:-1], gene_sequence[1:])])
        cpg_observed = gene_dinuc['CG']
        cpg_expected = sum([gene_dinuc[x] for x in cx_dinuc]) * sum([gene_dinuc[x] for x in gx_dinuc])
        cpg_expected = round(cpg_expected / len(gene_sequence), 2)
        cpg_bias = round(cpg_observed / cpg_expected, 4) if cpg_expected > 0 else 0
        print(gene, cpg_observed, cpg_expected, cpg_bias, sep='\t')
