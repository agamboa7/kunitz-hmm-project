# In Silico Modeling of the Kunitz-Type Domain with Profile Hidden Markov Model

This repository contains materials for the project developed for the *Laboratory of Bioinformatics I* course at the University of Bologna. The main objective is to build a Profile Hidden Markov Model (HMM) for the Kunitz-type protease inhibitor domain, starting from structural data and using it for functional annotation of protein sequences. The project involves retrieving data from the Protein Data Bank (PDB), sequence analysis, and the application of HMMs to search databases such as SwissProt. 
*Please note that this project is still in progress and the content in this repository is subject to change.

## Before Getting Started

Make sure you have [Conda](https://docs.conda.io/en/latest/) installed on your system. Then follow these steps to set up the environment and required tools:

```bash
# Create a Conda environment for the project
conda create --name lab_project

# Activate the environment
conda activate lab_project

# Install required bioinformatics packages
conda install -c bioconda cd-hit
conda install -c bioconda hmmer
conda install -c bioconda blast
conda install -c conda-forge biopython
```

## Step 1: Download PDB Files

Use the RCSB PDB Advanced Search with the following query:

**Query Conditions:**
- Data Collection Resolution ≤ 3.5 Å
- Pfam Annotation: `PF00014`
- Polymer Entity Sequence Length between 45 and 80 residues

Subsequently, generate a custom report with the following data:

**Custom Report Fields to Include:**
- **Identifier**: Entry ID  
- **Structure Data**: PDB ID, Data Collection Resolution  
- **Polymer Entity Data**: Sequence, Auth Asym ID, Annotation Identifier, Entity ID, Entry ID (Polymer Entity Identifiers)

> Export this report as a `.csv` file for the next step.


## Step 2: Convert Report to FASTA

Use the script `convert_to_fasta.sh` to convert the custom report into FASTA format.

```bash
bash convert_to_fasta.sh custom_report.csv
```
> Make sure to rename the PDBefold output file or name it appropriately so the script runs correctly.

This script takes the CSV file as an argument and generates a FASTA file called `kunitz_sequences.fasta` containing the protein sequences and identifiers. The script is included in this repository.

## Step 3: Cluster Sequences with CD-HIT

Run **CD-HIT** to cluster similar sequences using a 90% identity threshold (default). This helps reduce redundancy in your dataset.

Use the FASTA file generated in the previous step:

```bash
cd-hit -i kunitz_sequences.fasta -o kunitz_clustered.fasta -c 0.9
```
- `-i`: Input FASTA file  
- `-o`: Output file with clustered sequences  
- `-c`: Sequence identity threshold (`0.9` = 90%)

This step filters out highly similar sequences to retain representative diversity for profile HMM construction.

## Step 4: Extract Representative Sequences from Clusters

In the previous step, 160 sequences were processed with **CD-HIT**, resulting in **25 clusters**. Now, to retrieve the **representative sequences** (one per cluster), follow these steps:

1. **Extract the representative sequence IDs** from the `.clstr` file:
   ```bash
   grep '*' kunitz_clustered.fasta.clstr | cut -d '>' -f2 | cut -d '.' -f1 > representative_ids.txt
   ```
2. **Retrieve the corresponding sequences** from the original FASTA file and save them to a new file:
   ```bash
   for i in $(cat representative_ids.txt); do
         grep -A 1 "^>$i" kunitz_sequences.fasta | tail -n 2 >> kunitz_representatives.fasta
   done
   ```
3. **Filter Sequences by Length:** Filter out the short and long sequences from the representative FASTA file using the `count_and_filter.py` script provided in the repository. This script keeps only sequences between **55 and 85 amino acids** in length.

Run the script as follows:

```bash
python3 count_and_filter.py kunitz_representatives.fasta
```
After filtering, you should obtain a file `output_filtered.fasta` with 23 sequences that meet the length criteria.

## Step 5: Run Sequence Alignment on PDBefold

To analyze the filtered representative sequences in **PDBefold**, you need to provide a list of identifiers in the following format:
`PDBID:CHAIN`

To generate this list automatically from your filtered FASTA file, run the `pdb_efold_ids.sh` bash script included in the repository:

```bash
bash pdb_efold_ids.sh output_filtered.fasta
```

This will create a text file containing the correctly formatted IDs, ready for input into the PDBefold tool.

Subsequently, go to the [PDBefold website](https://www.ebi.ac.uk/msd-srv/ssm/) and fill out the **Multiple Alignment** submission form.

- Choose the option to **load a list of PDB codes**.
- Upload the file generated in the previous step (`pdb_efold_ids.txt`).
- Submit the form and wait a few moments for the tool to process the input.

Once the alignment is ready:

- **Download the FASTA file** provided by PDBefold.

The downloaded FASTA file may contain inconsistencies such as:
- Lowercase amino acids
- Extra information in the headers

To clean and standardize the file, use the `clean_fasta.sh` script provided in the repository:

```bash
bash clean_fasta.sh efold_output.txt
```
> Make sure to rename the PDBefold output file or name it appropriately so the script runs correctly.




