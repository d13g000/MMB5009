# MMB5009 Lectures

### Topic 1 - Sequencing and QC
#### **Exercise_1.2a**

This script trims DNA sequences and their corresponding quality scores in a FASTQ file. 
Users can specify whether to trim bases from the start, end, or both ends of each sequence.

<ins>Usage:<ins>

`python3 Exercise_1.2a.py <input_file> <mode> <n>`

    <input_file>: Path to the FASTQ file.
    <mode>: Trimming mode (start, end, or both).
    <n>: Number of bases to trim.

<ins>Example:<ins>

To trim 5 bases from the start of sequences in some_reads_1.fastq:

`python3 script1.py some_reads_1.fastq start 5`

<ins>Output:<ins>

The script creates a new FASTQ file named <input_file>_trimmed.fastq containing the trimmed sequences.


#### **Exercise_1.2b**

This script filters out sequences in a FASTQ file based on a user-provided mean quality score threshold.

<ins>Usage:<ins>

`python3 Exercise_1.2b.py <input_file> <quality_threshold>`

    <input_file>: Path to the FASTQ file.
    <quality_threshold>: The minimum mean quality score a sequence must have to be retained.

<ins>Example:<ins>

To filter some_reads_2.fastq and remove sequences with a mean quality below 30:

`python3 Exercise_1.2b.py some_reads_2.fastq 30`

<ins>Output:<ins>

The script filters out sequences with a mean quality score below the specified threshold.
It saves the filtered sequences in a new file named <input_file>_filtered.fastq.
It also displays statistics on the number of sequences removed and retained.



# MMB5009 Assignment 

## 1. Practical Task

### Filtering MS Data
#### clean_msdata.py

This script reads an Excel file containing mass spectrometry (MS) data and applies the following filters:
1. Maintains only rows where the Modifications column contains at least one 
   occurrence of \d+xMethyl, \d+xDimethyl, or \d+xTrimethyl 
   (case-insensitive), removing rows which contain any other type of 
   modification.
2. Maintains only rows where # Missed Cleavages > 0.
3. Excludes any peptides that have a modification on the first or last residue 
   (removes peptides with terminal PTMs).
4. Maintains only rows where the value in Abundance Ratio: (Cancer) / (Healthy) 
   is ≥ 1.2 (drops any row with missing or ≤ 1.2 ratio).

Furthermore, the script adds a new column named Abundance Ratio Flag (“High” 
or “Medium”) to all retained rows wherein:
    If the Abundance Ratio is ≥ 1.5 it is flagged as “High” and,
    If the Abundance Ratio is ≥ 1.2 and < 1.5 → it is flagged as “Medium”

The cleaned dataset is then written to a CSV.

<ins>Usage:<ins>

`python3 clean_msdata.py <input_excel> <output_csv>`

    <input_excel>: Path to the input Excel file (e.g. Healthy vs Cancer methyl peptide training set.xlsx).  
    <output_csv>: Path where the filtered CSV will be saved (e.g. cleaned_msdata.csv).  

<ins>Example:<ins>

To filter Healthy vs Cancer methyl peptide training set.xlsx and produce 
cleaned_msdata.csv:

`python3 clean_msdata.py Healthy vs Cancer methyl peptide training set.xlsx cleaned_msdata.csv`

<ins>Output:<ins>

A CSV file named cleaned_msdata.csv containing only those rows that meet all 
four criteria, plus a new column Abundance Ratio Flag wherein each retained 
row is labelled either "High" or "Medium" according to its abundance ratio.


### Annotate Peptides using UniProt
#### annotate_to_fasta.py

This script reads a filtered CSV (produced by clean_msdata.py) and, for each row:
1. Fetches the full UniProt sequence for the accession number in the Master 
   Protein Accessions column.
2. Parses the Positions in Master Proteins column to get the peptide’s 
   start and end locations within the master protein sequence taken from 
   UniProt.
3. Parses the Modifications in Master Proteins column to identify all 
   methyl‐type PTMs, distributing counts in such a way that:
      If the CSV file says 2xDimethyl [R14(100); R18(100)], ("Me2", 1) is assigned to positions 14 and 18.
      If the CSV file says 2xDimethyl [R585(100)], ("Me2", 2) is assigned to position 585.
      (i.e. the script is able to determine whether a modification occurs on two separate amino acids or on the same amino acid twice)
5. Builds a “pseudo-colored” protein sequence wherein:
      Residues outside the peptide are written in lowercase
      Residues inside the peptide are written in UPPERCASE
      If a residue inside is modified, <MeX:count> is appended to it immediately after that letter.
      Any PTMs falling outside the peptide window are ignored.
6. Writes one FASTA record per input row with the FASTA header including:
      UniProt accession number
      Peptide range
      Number of Methyl PTM sites within the length of the peptide
      Numeric abundance ratio (rounded to three decimals)
      Abundance Ratio Flag (High/Medium)
      Original UniProt FASTA header

<ins>Usage:<ins>

`python3 annotate_to_fasta.py <filtered_csv> <output_fasta>`

    <filtered_csv>: Path to the CSV from clean_msdata.py (e.g. cleaned_msdata.csv).  
    <output_fasta>: Path where the annotated FASTA will be written (e.g. annotated_cleaned_csv.fasta).

<ins>Example:<ins>

To annotate cleaned_msdata.csv and produce annotated_cleaned_csv.fasta:

`python3 annotate_to_fasta.py cleaned_msdata.csv annotated_cleaned_csv.fasta`

<ins>Output:<ins>

A FASTA file annotated_cleaned_csv.fasta where each entry looks like:

>P14678 | peptide=95-132 | mods_inside=2 | ratio=1.451 | flag=Medium | sp|P14678|...
mtvgksskmlqhidyrmrcilqdgrifigtfkafdkhmnlilcdcdefrkiVPLAGAAGGPGIGR<Me2:1>AAGRGIPAGVP... 

Where:
    The peptide region is shown in UPPERCASE, and each internal methylation site is tagged <MeX:count>,
    A blank line separates each FASTA record


### Generate Annotated CSV
#### generate_annotated_csv.py

This script reads the same filtered CSV (from clean_msdata.py) and produces 
a new CSV with the following columns - based on the researchers requirements:

1. Control/Cancer Sample Peptide: a multiline field showing /n
    `Control :`/n
    `<control peptide sequence>`/n
    
    `Cancer :`
    `<cancer peptide sequence>`
    control peptide sequence: extracted from Annotated Sequence column
    cancer peptide sequence: the same sequence with all methyl sites inserted as <MeX> (no count)
2. Modification (and which amino acid is methylated + Probability): the 
   exact methylation portion of the original Modifications column (e.g. 2xDimethyl [R14(100); R18(100)] or 1xTrimethyl [K585(100)]).
3. Positions in Master Proteins: copied verbatim from the input CSV and 
   representing the peptide's position within the master protein sequence.
4. Miscleavage: copied from # Missed Cleavages column.
5. Which Sample Replicate it was identified in: a multiline field showing the 
   sample IDs applied as:
    `Control:`
    `[S1]; [S2]; [S3];`

    `Cancer:`
    `[S4]; [S5]; [S6];`
    where the [S#] is only included if its corresponding Found in Sample: 
   [S#] column value ≠ “Not Found.” (if not found the sample number is not 
   included)
6. Abundance: a multiline field showing:
   `Control:`
   `<Abundances (Grouped): Healthy>`

    `Cancer:`
    `<Abundances (Grouped): Cancer>`
7. PSMs: copied from # PSMs column.

<ins>Usage:<ins>

`python3 generate_annotated_csv.py <filtered_csv> <output_csv>`

    <filtered_csv>: Path to the CSV from clean_msdata.py (e.g. cleaned_msdata.csv).  
    <output_csv>: Path where the new annotated CSV will be saved (e.g. annotated_msdata.csv).  

<ins>Example:<ins>

To create annotated_msdata.csv from cleaned_msdata.csv:

`python3 generate_annotated_csv.py cleaned_msdata.csv annotated_msdata.csv`

<ins>Output:<ins>

A CSV named annotated_msdata.csv with columns:
    Control/Cancer Sample Peptide,
    Modification (and which amino acid is methylated + Probability),
    Positions in Master Proteins,
    Miscleavage,
    Which Sample Replicate it was identified in,
    Abundance,
    PSMs

Where each row's cells contain the multiline fields as described above.
