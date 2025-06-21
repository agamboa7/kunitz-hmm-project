# In Silico Modeling of the Kunitz-Type Domain with Profile Hidden Markov Model

This repository contains materials for the project developed for the *Laboratory of Bioinformatics I* course at the University of Bologna. The main objective is to build a Profile Hidden Markov Model (HMM) for the Kunitz-type protease inhibitor domain, starting from structural data and using it for functional annotation of protein sequences. The project involves retrieving data from the Protein Data Bank (PDB) and UniProt, sequence analysis, and the evaluation of the performance of the HMM. 


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

## Step 5: Run Structure Alignment on PDBefold

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

## Step 6: Build the profile HMM from the cleaned alignment

Now that we have a cleaned multiple structure alignment, we can use it to build a profile HMM with `hmmbuild`. This model will represent the statistical profile of the Kunitz domain.

```bash
hmmbuild kunitz_domain.hmm aligned_clean.fasta
```
- `kunitz_domain.hmm`: Output file containing the HMM
- `aligned_clean.fasta`: Input multiple alignment (cleaned) in FASTA format

## Step 7: Prepare the datasets and create a BLAST database

To evaluate the performance of our HMM, we will need two datasets:

1. **All proteins *without* the Kunitz domain**  
2. **All proteins *with* the Kunitz domain**

You can retrieve these datasets directly from [UniProt](https://www.uniprot.org/) using the following advanced queries:

#### 🔹 Download datasets from UniProt:

- **All proteins *except* Kunitz domain:** NOT (xref:pfam-PF00014) AND (reviewed:true)
- **All proteins *with* Kunitz domain:** (xref:pfam-PF00014) AND (reviewed:true)

Download the results in **FASTA format**.

To prepare for filtering and classification, extract the UniProt accession IDs from the FASTA headers of both the Kunitz and non-Kunitz datasets:

```bash
# Extract IDs from reviewed Kunitz domain sequences
grep ">" uniprotkb_xref_pfam_PF00014_AND_reviewe_2025_06_09.fasta | cut -d "|" -f2 > all_kunitz.ids
```

```bash
# Extract IDs from reviewed non-Kunitz sequences
grep ">" uniprotkb_NOT_xref_pfam_PF00014_AND_rev_2025_06_09.fasta | cut -d "|" -f2 > all_not_kunitz.ids
```
> Make sure to use the adequate names of both UniProt files to ensure the scripts run correctly.

#### 🔹 Create a BLAST database from the Kunitz protein sequences:

Once you've downloaded the Kunitz protein sequences (with Pfam: PF00014), use the following command to create a BLAST database:

```bash
makeblastdb -in uniprotkb_xref_pfam_PF00014_AND_reviewe_2025_06_09.fasta -dbtype prot -out kunitz_db
```
- `-in`: Input FASTA file containing Kunitz domain proteins
- `-dbtype prot`: Specifies that the database is protein-based
- `-out`: Name of the output BLAST database

## Step 8: Filter similar sequences to avoid bias

To ensure that our evaluation isn't biased by sequences that are too similar to those used to build the HMM, we will perform sequence similarity filtering.

#### 1. Run BLAST between training sequences and UniProt Kunitz set

Use the clean aligned FASTA (from Step 6) as the query against the Kunitz domain BLAST database:

```bash
blastp -query aligned_clean.fasta -db kunitz_db -out kunitz_filter.blast -outfmt 7
```
- `-query`: PDB-based HMM training sequences
- `-db`: Kunitz BLAST database created in Step 7
- `-out`: Output file with the filtering
- `-outfmt`: Tabular output with comment lines (easier for filtering)

#### 2. Identify and filter highly similar sequences

To remove sequences that are highly similar (≥95% identity and alignment length ≥50 residues), run:

```bash
grep -v "^#" kunitz_filter.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u > to_remove.txt
```

#### 3. Clean the list to keep only UniProt IDs

Remove extra characters and keep just the UniProt IDs with the following command:

```bash
cut -d '|' -f2 to_remove.txt > to_remove.ids
```

#### 4. Filter out those IDs from the full Kunitz set

Filter these IDs from the all_kunitz.ids file, which contains the original UniProt IDs from the full Kunitz set: 

```bash
comm -23 <(sort all_kunitz.ids) <(sort to_remove.ids) > cleaned_kunitz_ids.txt
```

> This step reduces redundancy: The total number of IDs went from 398 original Kunitz IDs to 365 non-redundant entries.

## Step 9: Generate random positive and negative datasets for HMM evaluation

To evaluate the HMM classifier fairly, generate two random subsets from the positive (Kunitz) and negative (non-Kunitz) datasets:

```bash
# Shuffle the cleaned Kunitz (positive) IDs
sort -R cleaned_kunitz_ids.txt > cleaned_kunitz_random.txt
```

```bash
# Split into two positive sets
head -n 199 cleaned_kunitz_random.txt > pos_1.txt
tail -n 199 cleaned_kunitz_random.txt > pos_2.txt
```

```bash
# Shuffle the non-Kunitz (negative) IDs
sort -R all_not_kunitz.ids > all_not_kunitz_random.txt
```

```bash
# Split into two negative sets
head -n 286416 all_not_kunitz_random.txt > neg_1.txt
tail -n 286416 all_not_kunitz_random.txt > neg_2.txt
```

## Step 10: Extract FASTA sequences for each set using `get_seq.py`

The `get_seq.py` script (included in the repository) retrieves FASTA sequences by matching IDs from a text file against a larger FASTA file. Use this script to extract the sequences for each of your datasets:

```bash
# Extract positive sequences from Kunitz FASTA
python3 get_seq.py pos_1.txt uniprotkb_xref_pfam_PF00014_AND_reviewe_2025_06_09.fasta
python3 get_seq.py pos_2.txt uniprotkb_xref_pfam_PF00014_AND_reviewe_2025_06_09.fasta
```

```bash
# Extract negative sequences from non-Kunitz FASTA
python3 get_seq.py neg_1.txt uniprotkb_NOT_xref_pfam_PF00014_AND_rev_2025_06_09.fasta
python3 get_seq.py neg_2.txt uniprotkb_NOT_xref_pfam_PF00014_AND_rev_2025_06_09.fasta
```

> Each command will produce a FASTA file with the same base name as the ID file, but with a `_f.fasta` extension (e.g., `pos_1_f.fasta`, `neg_1_f.fasta`, etc.). These files will be used in the next steps to classify and evaluate sequences with the HMM.

## Step 11: Search sequences with your HMM and prepare classification files

Use the `kunitz_domain.hmm` HMM profile to scan each of the data folds using `hmmsearch`. The `--tblout` option saves a summary table of the hits. The `-Z 1000` normalizes the database size, and `--max` reports the best domain only.

```bash
# Run hmmsearch on each set
hmmsearch -Z 1000 --max --tblout pos_1.out kunitz_domain.hmm pos_1_f.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out kunitz_domain.hmm pos_2_f.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out kunitz_domain.hmm neg_1_f.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out kunitz_domain.hmm neg_2_f.fasta
```

Now convert the output files to `.class` format for performance evaluation:

```bash
# Transform hmmsearch output to .class format
grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > neg_2.class
```

Finally, merge the positive and negative sets into two main classification files:

```bash
cat pos_1.class neg_1.class > set_1.class
cat pos_2.class neg_2.class > set_2.class
```

## Step 12: Evaluate the model performance with different E-value thresholds

Now that both `set_1.class` and `set_2.class` have been created, evaluate the performance of the HMM using the `performance.py` script (available in the repository). This script takes a `.class` file and an E-value threshold as input.

In this step, we iterate through E-value thresholds from `1e-1` to `1e-12` (covering a range of stringency levels) to assess how strict or relaxed thresholds affect classification performance.

Run the script for **set 1**:

```bash
for i in $(seq 1 12); do
    echo "Threshold: 1e-$i"
    python3 performance.py set_1.class 1e-$i
    echo ""
done
```

Now repeat the same process for **set 2**:

```bash
for i in $(seq 1 12); do
    echo "Threshold: 1e-$i"
    python3 performance.py set_2.class 1e-$i
    echo ""
done
```

For each threshold (e.g., 1e-3, 1e-4, etc.), the `performance.py` script:

- Builds a confusion matrix.
- Calculates performance metrics such as:
  - `q2`: overall accuracy
  - `MCC`: Matthews Correlation Coefficient
  - `TPR`: true positive rate (sensitivity/recall)
  - `PPV`: positive predictive value (precision)

It evaluates the classification performance **twice**:
- Once using the **full-sequence E-value** (column 3 in the `.class` file).
- Once using the **best domain E-value** (column 4 in the `.class` file).

This dual evaluation allows to compare whether using the full-sequence E-value or the best domain hit is more accurate for classifying.

📌 Note: A third argument can also be specified when runnin the `performance.py` script:
- `1` → only evaluate **full sequence**
- `2` → only evaluate **best domain**
- `0` (or omitted) → evaluate **both**

## Step 13 (Optional): Visualize the results using R

To visualize how well the HMM performs at different thresholds, we can generate a plot of **MCC vs. E-value threshold** for both datasets and evaluation types (full sequence vs. best domain).

#### 1. Generate a `.tsv` file with performance results

Instead of printing the metrics to the terminal, run the `performance_tsv.py` script included in this repository. It works the same way as `performance.py`, but stores the output in a structured `.tsv` format, and does the threshold iteration inside the python script.

Run it like this:

```bash
python3 performance_tsv.py set_1.class
python3 performance_tsv.py set_2.class
```
#### 2. Visualize the results in RStudio

After obtaining the TSV files, run the `performance_graph.R` script (included in this repository) in RStudio to create a plot comparing the Matthews Correlation Coefficient (MCC) values for both sets.

The script used in our project loads both TSV files and plots them together with:
- Set 1 shown in purple
- Set 2 shown in hotpink
- Line types for Full Sequence vs. Best Domain

This allows to visually assess:
- How MCC changes with stricter or looser E-value thresholds.
- Whether using full sequence or best domain produces more stable and accurate classification.

> In this case (first run for Laboratory of Bioinformatics I project), all MCC values were high (above 0.92), with Full Sequence reaching above 0.99 but also showing slightly more variability than Best Domain. Results and Discussion will be developed in the Project Report.


✅ This concludes the pipeline!

---

### Questions or Feedback?

If you have any questions, encounter issues, or would like to contribute to improving this pipeline, feel free to contact:

#### Andrea Gamboa  |  **Email:** andrea.arriolagamboa@studio.unibo.it
Department of Pharmacy and Biotechnology Alma Mater Studiorum – Università di Bologna 

> This project was built for academic purposes as part of the Laboratory of Bioinformatics I course at the University of Bologna (A.Y. 2024/2025)

