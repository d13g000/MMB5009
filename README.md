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

<ins>Dependencies:<ins>

- Python 3+
- pandas
- requests (for annotate_to_fasta.py)
- openpyxl (for clean_msdata.py to read Excel)

Install via:

`conda install pandas requests openpyxl`

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
1. Control/Cancer Sample Peptide: a multiline field showing 

    `Control :`

    `<control peptide sequence>`
    
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


### Example Workflow

1. Filter raw MS data 

`python3 clean_msdata.py Healthy vs Cancer methyl peptide training set.xlsx cleaned_msdata.csv`

2. Build annotated FASTA

`python3 annotate_to_fasta.py cleaned_msdata.csv annotated_cleaned_csv.fasta`

3. Generate final annotated CSV

`python3 generate_annotated_csv.py cleaned_msdata.csv annotated_msdata.csv`



## 2. Generating a Dataset of Protein Data

<ins>Dependencies:<ins>

- Python 3+
- pandas
- requests
- openpyxl

Install via:

`conda install pandas requests openpyxl`

### Generating Protein Dataset
#### generate_protein_dataset.py

This script reads a CSV or Excel file containing UniProt accessions in a 
specified column, deduplictes them (identifies and removes duplicate protein 
accessions), prompts the user to process either all or the first N 
accessions, fetches detailed annotations via the EBI Proteins 
API, and outputs a structured JSON file.
1. Opens and reads the CSV/Excel file, pulling values from the --column 
   header (defaults to "Master Protein Accessions"). In doing so the script 
   splits cells on semicolons (;), trims the whitespace, and builds a sorted 
   list of unique accession numbers. 
2. Lists the number of unique protein accessions found and estimates the 
   runtime to process them all (number of accessions x delay time set by 
   user), then allows the user to choose between processing: all, the first N 
   (where N is entered by the user), or aborting.
3. If the user requests processing the script fetches the protein accession 
   number using https://www.ebi.ac.uk/proteins/api/proteins/{accession} and 
   parses out the respective protein's:
        - Existence and Version
        - Protein Name
        - Accessions (primary/secondary)
        - Gene Name and Symbol
        - Status (Functions, Interactions, etc.)
        - Organism (Taxonomy number, Scientific name, Common name)
        - Variants and Variant IDs
        - Regions (all non-variant features)
        - GO Annotations
        - Associated Diseases (Texts, Evidences, Scope)
        - Families and Domains
        - Sequence (Canonical sequence)
        - Isoforms (IDs, Name, Status, Sequence)
4. Writes a JSON file.

<ins>Usage:<ins>

`python3 generate_protein_dataset.py \
  --input <input_csv> \
  --column "Master Protein Accessions" \
  --output <output_json> \
  [--delay 0.1]`

    --input <input_csv>: Path to the CSV/Excel containing column with 
    accession numbers (e.g. Healthy vs Cancer methyl peptide training set.csv).
    --column: Header name (defaults to "Master Protein Accessions".
    --output <output_json>: Path where the JSON file will be saved 
    (e.g. protein_database.json).
    --delay: Time (in seconds) between API calls (defaults to 0.1).

<ins>Example:<ins>

To create protein_dataset.json from Healthy vs Cancer methyl peptide 
training set.csv:

`python3 generate_protein_dataset.py \
    --input Healthy vs Cancer methyl peptide training set.csv \
    --column "Master Protein Accessions" \
    --output protein_dataset.json \
    --delay 0.1`

<ins>Output:<ins>

A JSON file containing a list of proteins and their respective parsed data.


## 3. Implementing the Bad Character and Good Suffix Rules for Exact Matching Algorithms

<ins>Dependencies:<ins>

- Python 3+
- other libraries used are in the Python standard library

### Implementing an Exact Matching Algorithm (Boyer-Moore using both Bad Character and Good Suffix Rules, including the option of the Pigeon-Hole Principle)
#### bm_search.py

This script loads a FASTA reference sequence and implements exact and 
approximate matching algorithms. Furthermore, it also benchmarks all methods 
and reports comparative information.
1. Exact Matching:
   This script runs exact matches using two separate algorithms:
      - Naive sliding window search (one amino acid at a time).
      - Boyer-Moore with both Bad Character (skip to next matching character 
        when encountering a mismatch) and Good Suffix (skip to next matching 
        motif when encountering a mismatch) heuristics. In doing so the script 
        chooses the method (Bad Character/Good Suffix) that gives the largest 
        shift on mismatch, ensuring that it handles overlaps within the 
        sequence correctly.
2. Approximate Matching (<= N mismatches) via the Pigeon-Hole Principle:
   - During this algorithm, the pattern input by the user is (1) initially 
     split into N+1 pieces, (2) next an exact search of each split piece is 
     conducted using the Boyer-Moore algorithm, and (3) finally candidate 
     alignments are verified based on Hamming distance.
3. Benchmarking:
   - While benchmarking the script shows: the number of hits found by each 
     of the three algorithms (Naive/Boyer-Moore/Pigeon-Hole), the number of 
     alignments inspected (showing how many the Boyer-Moore method skips 
     when compared to the Naive one), and the elapsed time (formatted into 
     minutes and seconds if more than 60s).

<ins>Usage:<ins>

`python3 bm_search.py <reference_fasta> <pattern> [--mismatches N] 
[--benchmark]`

    <reference_fasta>: Path to the FASTA file with the query sequence (e.g. 
    chr1.fa).
    <pattern>: Sequence to search for.
    --mismatches N: The number of mismatches allowed via the Pigeon-Hole 
    principle (defaults to 0).
    --benchmark: Report method-specific information.

<ins>Example:<ins>

To run Boyer-Moore and Pigeon-Hole algorithms to match a user-specified DNA 
sequence within Chromosome 1: 

`python3 bm_search.py chr1.fa AACCGGAATTACGTA --mismatches 2 --benchmark`

<ins>Output:<ins>

Terminal output outlining the comparison of each method:
- Number of hits and alignments inspected by each method,
- Number of alignments skipped using Boyer-Moore when compared to the Naive 
  method,
- The elapsed time for each method (in human-readable format if above 60s),
- A list of detailed match positions including the number of mismatches.


    Naive (exact):          {X} hits, inspected {Y} alignments, time = {Z}
    Boyer-Moore (exact):    {X} hits, inspected {Y} alignemnts, time = {Z}
    
    BM skipped {N} alignments compared to naive

    Pigeon-Hole (≤{n} mm):  {X} hits, time = {Z}

    Found {X} hits (pos -> mismatches):
    {list of positions} {list of mismatches}


## 4. Mini BLAST

<ins>Dependencies:<ins>

- Python 3+
- biopython

Install via:

`conda install biopython`

### Implementing a Simplified Version of BLAST
#### mini_blast.py

This script implements a simplified, configurable BLAST-style search to rank 
a small database of protein sequences by similarity to a query sequence. The 
user is able to choose between:
- Ungapped Highest-Scoring Pair (HSP) Extension using k-mer word size seeding, 
exact matches, and X-drop ungapped extension with BLOSUM62 scoring.
- Gapped local alignment (Smith-Waterman style) using Biopython's pairwise2 
with BLOSUM62 and affine gap penalties.
After running, the script returns information from the query's FASTA header, 
the mode of alignment used, and the ranking of target sequences ranging from 
most similar to least similar based on their HSPs or local alignment scores 
(depending on alignment method chosen).

<ins>Features:<ins>

1. Configurable word size (-k / --wordsize) for seeding when choosing the 
   ungapped mode of alignment (defaults to 3).
2. Choice of ungapped (--ungapped) and gapped (--gapped) alignment modes 
   (defaults to ungapped).
3. BLOSUM62 substitution matrix for scoring matches/mismatches.
4. Affine gap penalties when using gapped mode where:
   - gap open = -10, 
   - gap extend = -1.

<ins>Usage:<ins>

`python3 mini_blast.py [--wordsize N] [--gapped | --ungapped] <query.fa> 
<target.fa>`

    --wordsize: The word size for k seeding (only used in ungapped mode, 
    defaults to 3).
    --gapped: The option to use ungapped HSP extension (default mode).
    --ungapped: The option to use gapped local alignment (Smith-Waterman 
    method).
    <query.fa>: FASTA file containing single query sequence.
    <target.fa>: FASTA file containing database of sequences to compare to 
    query sequence.

<ins>Example:<ins>

To run ungapped alignment with default k-mer word size seeding:

`python3 mini_blast.py query.fa target.fa`

To run ungapped alignment with k-mer word size = 5:

`python3 mini_blast.py -k 5 --ungapped query.fa target.fa`

To run gapped local alignment (Smith-Waterman method):

`python3 mini_blast.py --gapped query.fa target.fa`

<ins>Output:<ins>

Terminal output outlining the query sequence and its associated description, 
the mode of alignment used, and the ranking of the target sequences based on 
alignment score.


    Query: {query sequence ID}
    Description: {query sequence description}
    Mode: {alignment algortihm used} {word size seed [if using --ungapped]}

    Rank        Target_ID       Score       Description
    {ranked list of database sequences and associated information}


## 5. RNA-Seq Analysis Pipeline - Normalisation

<ins>Dependencies:<ins>

- Python 3+
- pysam

Install via:

`pip install pysam`

### Normalisation during RNA-Seq Analysis
#### normalise_rpm_rpkm.py

This command-line utility computes Reads Per Million (RPM) and Reads Per 
Kilobase Million (RPKM) from gene counts generated by featureCounts output 
(from the subread software). It parses the featureCounts file for gene 
lengths and raw counts, tallies total mapped reads in a coordinate-sorted BAM, 
and outputs a TSV table with normalisation metrics.
1. Total Mapped Reads: Firstly, the script excludes unmapped, secondary, 
   QC-failed, and duplicate reads using pysam.
2. featureCounts Parsing: Secondly, the script reads and parses the Geneid, 
   Length, and the last column as raw counts (sample_chr21_22_SE_mapped.bam).
   

    Geneid   Chr   Start   End   Strand   Length   sample_chr21_22_SE_mapped.bam (raw counts information)

3. Normalisation Metrics: Next, the script runs RPM and RPKM normalisations;
- RPM: Raw counts per million total reads,
- RPKM: RPM per kilobase of region length.
5. Progress Reporting and Robust Error Handling: Finally, the script 
   outputs its progress and skips malformed lines, outputting any errors 
   encountered.

<ins>Usage:<ins>

`python3 normalise_rpm_rpkm.py \
    --counts path/to/counts/file.txt \
    --bam    path/to/alignments/file.bam \
    --out    path/to/output/file.tsv`

    -c, --counts: featureCounts output TSV file (having format shown in 2.).
    -m, --bam: Mapped BAM file (coordinate-sorted and indexed).
    -o, --out: Output TSV file.

<ins>Example:<ins>

`python3 normalise_rpm_rpkm.py \
    --counts gene_counts_chr21_22.txt \
    --bam    sample_chr21_22_SE_mapped.bam \
    --out    gene_norm_chr21_22.tsv`

<ins>Output:<ins>

Running the script returns terminal output showing its progress:

    Counting total mapped reads...
    Total mapped reads: 11407404
    Parsing featureCounts output...
    Computing RPM/RPKM and writing output...
    Done. Results in gene_norm_chr21_22.tsv

Upon completion, the TSV file generated outlines the gene's ID, length, 
number of raw counts, and normalised counts generated by RPM and RPKM 
normalisations.

    gene_id     length_bp       raw_counts      RPM                 RPKM
    {geneID}    {length}        {counts}        {RPM_norm_value}    {RPKM_norm_value}