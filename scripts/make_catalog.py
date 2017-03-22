#! /usr/bin/env python

import pandas
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_catalog', required=True, type=str, dest='incat',
                    help='The input catalog (e.g., COSMOS.in).')

parser.add_argument('--z_catalog', required=True, type=str, dest='zcat',
                    help='The redshift catalog, output from photo-z.')

parser.add_argument('--out_catalog', required=True, type=str, dest='outcat',
                    help='The final, combined output catalog.')

args = parser.parse_args()

zcat = pandas.read_table(args.zcat, delim_whitespace=True)

incat = pandas.read_csv(args.incat, delim_whitespace=True,
                        names=["id",
                               'ra',
                               'dec',
                               'Imag',
                               'zin',
                               'template',
                               'ebv',
                               'rlaw',
                               'scale',
                               'line-scale'])

outcat = pandas.concat([incat, zcat], axis=1)

outcat.to_csv(args.outcat, sep=' ', index=False)

