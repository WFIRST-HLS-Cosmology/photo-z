#! /usr/bin/python

import sys

if len(sys.argv) < 3 or len(sys.argv) > 4:
    sys.stderr.write("\nUsage: " + sys.argv[0] + " column_number input_file [output_file]\n\n")
    sys.exit(1)

col = int(sys.argv[1]) - 1 # the first column is column 1---not column 0

if col < 0:
    sys.stderr.write("Error: The column number must be 1 or greater\n")
    sys.exit(1)

filename = sys.argv[2]

input_file = open(filename)

contents = []

for line in input_file:
    contents.append(line.split())

if col >= len(contents[0]):
    sys.stderr.write("Error: The column number is out of range (more than the number of columns in the input file).\n")
    sys.exit(1)

is_int = False
is_float = False

try:
    x = int(contents[0][col])
    is_int = True
except ValueError:
    is_int = False

try:
    x = float(contents[0][col])
    is_float = True
except ValueError:
    is_float = False

if is_float:
    contents.sort(key=lambda entry : float(entry[col]))
elif is_int:
    contents.sort(key=lambda entry : int(entry[col]))
else:
    contents.sort(key=lambda entry : entry[col])

if len(sys.argv) == 3:
    for row in contents:
        line = " ".join(row)
        sys.stdout.write(line + "\n")
else:
    output_file = open(sys.argv[3], 'w')
    for row in contents:
        line = " ".join(row)
        output_file.write(line + "\n")