#!/bin/bash

###############################################################################
# Script Name: clean_fasta.sh
# Description: Takes a FASTA alignment file with headers in the format
#              >PDB:CHAIN + info and produces a cleaned version where:
#              - IDs are reformatted to >PDB_CHAIN (uppercase)
#              - Sequences are fully capitalized
# Usage:       ./clean_fasta.sh input.fasta
# Output:      aligned_clean.fasta
# Author:      Andrea Arriola Gamboa
# Date:        26-05-2025
###############################################################################

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "ERROR:"
    echo -e "This program takes one argument from the command line: input file name.\nExiting..."
    exit 1
fi

input_file="$1"
output_file="aligned_clean.fasta"

# Process the FASTA file
awk '/^>/ {
    split($1, id, ":");
    print ">" toupper(id[2]) "_" toupper(id[3])
}
!/^>/ {
    print toupper($0)
}' "$input_file" > "$output_file"

echo "Cleaned FASTA file saved as: $output_file"
