#!/bin/bash
#
# This script takes an input file and dump it in hexadecimal (1 byte) format
# The second input is the name of the output file
#
# The produced file can be included in C++ source code as an array of char:
#
#   char inputString[] = {
#       #include "<outputfile>"
#       0x00
#   };

inputfile="$1"
outputfile="$2"

# 1. gzip output file
# 2. convert to hex
# 3. drop last line of output
# 4. print all but first parameter (offset) in C++ formatting

cat "$inputfile" | \
gzip -9 | \
od -v -t x1 | \
awk 'NR>1{print buf}{buf = $0}' | \
awk '{for (i=2; i<=NF; i++) print "0x"$i","}' > "$outputfile"
