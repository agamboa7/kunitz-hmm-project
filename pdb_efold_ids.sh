#!/bin/bash

###############################################################################
# Script Name: pdb_efold_ids.sh
# Description: Extracts PDB IDs with chain identifiers from a FASTA file
#              (formatted as >PDBID_CHAIN) and outputs them in the format
#              PDB:CHAIN (e.g., 1ABC:A) to a text file.
# Usage:       ./pdb_efold_ids.sh input.fasta
# Output:      pdb_efold_ids.txt
# Author:      Andrea Arriola Gamboa
# Date:        26-05-2025
###############################################################################

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "ERROR:"
    echo -e "This program takes one argument from the command line: input file name.\nExiting..."
    exit 1
fi

input_file="$1"
output_file="pdb_efold_ids.txt"

# Extract PDB IDs with chain identifiers in the format PDB:CHAIN
awk 'substr($0,1,1) == ">" {
    id = substr($0, 2)
    split(id, a, "_")
    print a[1] ":" a[2]
}' "$input_file" > "$output_file"

echo "IDs written to $output_file"
