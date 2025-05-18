#!/usr/bin/env python3

"""
> natural_sort.py <

Does it says on the tin: can be imported into other scripts for the
natural_sort function, or can be called from the command line to sort an input
list.

Natural sort: 15, 9, 101 --> 9, 15, 101
UNIX sort:    15, 9, 101 --> 101, 15, 9
"""
import argparse
import re
import sys

def natural_sort(input_list, reversed=False):
    """Sort a list using natural sort order."""
    def chunked_text(x):
        return [int(y) if y.isdigit() else y for y in re.split('([0-9]+)', x)]
    return sorted(input_list, key=chunked_text, reverse=reversed)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Script takes in an input file, natural sorts it, and spits it back out.""")
    parser.add_argument('unsorted', metavar='text_file',
                        type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin, help='unsorted file.')
    parser.add_argument('-r', '--reverse', action='store_true', default=False,
                        help='reverse the sort order.')

    args = parser.parse_args()

    # Read, sort, and print lines
    sorted_lines = natural_sort((line.strip() for line in args.unsorted), args.reverse)
    sys.stdout.write("\n".join(sorted_lines) + "\n")
