# MB5009

**Exercise_1.2a**

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


**Exercise_1.2b**

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