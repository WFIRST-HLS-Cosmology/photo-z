#! /usr/bin/env python

"""
filter an input stream by only printing lines in which the value of a 
specified column is less than a specified threshold value.

This has the same effect as the awk filter:

   awk '{if ($column_index < threshold) print $0}'

"""
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--column', required=True, type=int, dest='column_idx',
                    help='The column number of interest, using one-based indexing.')

parser.add_argument('--limit', required=True, type=float, dest='max_val',
                    help='The maximum value to allow through the filter.')

args = parser.parse_args()


for line in sys.stdin:
    column_value = float(line.split()[args.column_idx - 1])
    if column_value <= args.max_val:
        sys.stdout.write(line)
