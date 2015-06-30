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

cat "$inputfile" | \
gzip -9 | \
od -v -t x1 | \
sed \
's/^[0-9A-Fa-f]\{1,\} *//
/^ *$/d
s/ \{1,\}/,0x/g
s/^/0x/
s/$/,/' > "$outputfile"
