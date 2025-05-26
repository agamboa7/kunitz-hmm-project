#!/usr/bin/env python3
# Author: Andrea Arriola Gamboa
# Python version: 3.12.9

"""
count_and_filter.py

This script prints the length of each (single-line) sequence in a FASTA file.
It then prompts the user to confirm whether they want to filter the sequences
by length (between 55-85). If the user types 'yes', an output file is generated
containing only the filtered sequences.

Usage:
    python count_and_filter.py <input_file.fasta>
"""

from sys import argv

def count_aa_in_fasta(file):
    """
    Reads a FASTA file and prints the length of each
    sequence with its corresponding identifier.
    Then it prints the number of sequences in the file.

    Parameters:
    file (str): FASTA file to be processed.
    """
    with open(file, "r") as fasta_file:
        current_id = ""
        num_sequences = 0
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                current_id = line[1:]
                num_sequences = num_sequences + 1
            else:
                aa_count = len(line)
                print(f"{current_id} -> {aa_count}")
        print(f"\nNumber of sequences in the input file: {num_sequences}")
            

def filter_fasta_by_len(input_file, output_file, min_len=55, max_len=85):
    """
    Filters sequences from a FASTA file by length and writes the valid
    sequences to a new output file.
    Then it prints the number of sequences in the output file.

    Parameters:
    input_file (str): Input FASTA file.
    output_path (str): Output file for filtered sequences.
    min_len (int): Minimum allowed sequence length.
    max_len (int): Maximum allowed sequence length.
    """
    with open(input_file, "r") as input_file, open(output_file, "w") as output_file:
        current_id = ""
        current_seq = ""
        for line in input_file:
            line = line.rstrip()
            if line.startswith(">"):
                # Check and write the previous sequence if within range
                if current_id and min_len <= len(current_seq) <= max_len:
                    output_file.write(f">{current_id}\n")
                    output_file.write(f"{current_seq}\n")
                current_id = line[1:]
                current_seq = ""
            else:
                current_seq += line
        # Check the last sequence
        if current_id and min_len <= len(current_seq) <= max_len:
            output_file.write(f">{current_id}\n")
            output_file.write(f"{current_seq}\n")

# Run the program
if __name__ == "__main__":
    if len(argv) != 2:
        print("ERROR: This program takes a FASTA file as an argument.")
        print("Usage: python seq_checker.py <input_file.fasta>")
        exit()
    else:
        input_file = argv[1]
        output_file = "output_filtered.fasta"

        # Step 1: Print all sequence lengths
        print("Sequence lengths:")
        count_aa_in_fasta(input_file)

        # Step 2: Ask user if they want to filter
        user_choice = input("Do you want to filter sequences with length 55–85? (yes/no): ").strip().lower()
        if user_choice == "yes":
            filter_fasta_by_len(input_file, output_file, min_len=55, max_len=85)
            print(f"\nFiltered sequences written to: {output_file}")
        # If they filter, the program also prints the new number of sequences
            with open(output_file,"r") as output_file:
                num_sequences = 0
                for seq in output_file:
                    if seq.startswith(">"):
                        num_sequences = num_sequences + 1
                print(f"\nNumber of sequences in the output file: {num_sequences}")
        elif user_choice == "no":
            print("\nNo filtering performed.")
        else:
            print("\nInvalid option. No filtering performed.")
