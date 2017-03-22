#! /usr/bin/env python

import sys

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " sorted_cosmos_catalog"
    exit(1)

catalog_filename = sys.argv[1]

catalog = open(catalog_filename)

redshift_threshold = 0.0005

part = 0

n_redshifts = 0

max_redshifts = 20

prev_z = float(catalog.readline().split()[6])

catalog.seek(0)

current_file = open(catalog_filename + "-part-" + str(part), 'w')

for line in catalog:
    redshift = float(line.split()[6])

    if abs(redshift - prev_z) > redshift_threshold:
        n_redshifts += 1
        prev_z = redshift

    if n_redshifts == max_redshifts:
        current_file.close()
        part += 1
        n_redshifts = 0
        current_file = open(catalog_filename + "-part-" + str(part), 'w')

    current_file.write(line)

current_file.close()