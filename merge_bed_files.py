#!/usr/bin/env python3

"""
> merge_bed_files.py <

Given modkit traditional bed files, combine the meth and valid reads
together, compute the meth %.
"""
import argparse
import csv
import sys
from collections import defaultdict

import natural_sort

parser = argparse.ArgumentParser(description="""
Given filtered modkit traditional bed files, combine the meth and 
valid reads together, compute the meth %.""")

parser.add_argument('bed_files', metavar="bed_files",
                    type=argparse.FileType('r'), nargs='+',
                    help="List of bed files from modkit.")
parser.add_argument('-v', action='store_true',
                    help="verbose mode, prints progress to stderr.")
args = parser.parse_args()

# Use defaultdict for nested dictionaries
combined_data = defaultdict(lambda: defaultdict(lambda: [0, 0]))
counter_rows = 0

for c in args.bed_files:
    tsv_reader = csv.reader(c, delimiter='\t')
    for row in tsv_reader:
        if not row:
            continue

        if args.v:
            counter_rows += 1
            if counter_rows % 1000000 == 0:
                print(f'{counter_rows} rows processed...', file=sys.stderr)

        scaf = row[0]
        pos = int(row[1])
        meth = int(row[11])
        valid = int(row[9])

        # Update meth and unmeth counts
        combined_data[scaf][pos][0] += meth
        combined_data[scaf][pos][1] += valid

# Prepare output
output_writer = csv.writer(sys.stdout, delimiter='\t')
for scaf in natural_sort.natural_sort(combined_data.keys()):
    for pos in sorted(combined_data[scaf]):
        meth, valid = combined_data[scaf][pos]
        meth_pct = round(meth / valid * 100, 4)
        output_writer.writerow([scaf, pos, pos, meth_pct, meth, valid])
