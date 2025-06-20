"""
normalize_rpm_rpkm.py

Compute RPM and RPKM values from featureCounts output.
Follows a structured pipeline:
  1) Parse command-line arguments and options
  2) Count total mapped reads from BAM (excluding unmapped, secondary,
  qcfail, duplicates)
  3) Load featureCounts output to extract gene lengths and raw counts
  4) Compute RPM (reads per million) and RPKM (reads per kilobase per million)
  5) Write normalized values to output TSV

Assumptions and options:
  - Input featureCounts TSV contains header with 'Geneid' and 'Length' columns,
    and exactly one sample count column at the end.
  - BAM file must be coordinate-sorted and indexed.
  - RPM and RPKM are reported with six decimal places.
"""

import argparse
import sys
import pysam

# -------- Step 1: Count total mapped reads ----------
def count_mapped_reads(bam_path):
    """
    Counts total mapped reads in a BAM file, excluding unmapped, secondary,
    QC failed, and duplicate reads.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    total = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate:
            continue
        total += 1
    bam.close()
    return total

# -------- Step 2: Load featureCounts output ----------
def load_featurecounts(counts_path):
    """
    Parses featureCounts output, extracting gene lengths and raw counts.
    """
    gene_lengths = {}
    raw_counts = {}
    with open(counts_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "Geneid":
                idx_len = parts.index("Length")
                idx_cnt = len(parts) - 1
                continue
            gene = parts[0]
            length = int(parts[idx_len])
            count = int(parts[idx_cnt])
            gene_lengths[gene] = length
            raw_counts[gene] = count
    return gene_lengths, raw_counts

# -------- Step 3: Compute RPM and RPKM ----------
def compute_normalization(gene_lengths, raw_counts, total_mapped):
    """
    Computes RPM and RPKM values for each gene.
    """
    million_factor = total_mapped / 1e6 if total_mapped > 0 else 1
    results = []
    for gene, length in gene_lengths.items():
        count = raw_counts.get(gene, 0)
        rpm = count / million_factor if million_factor > 0 else 0
        rpkm = rpm / (length / 1000) if length > 0 else 0
        results.append((gene, length, count, rpm, rpkm))
    return results

# -------- Step 4: Write output TSV ----------
def write_output(results, out_path):
    """
    Writes normalised counts (RPM, RPKM) to a TSV file.
    """
    with open(out_path, 'w') as out:
        out.write("gene_id\tlength_bp\traw_counts\tRPM\tRPKM\n")
        for gene, length, count, rpm, rpkm in results:
            out.write(f"{gene}\t{length}\t{count}\t{rpm:.6f}\t{rpkm:.6f}\n")

# -------- Main entry point ---------
def main():
    parser = argparse.ArgumentParser(
        description="Compute RPM and RPKM from featureCounts output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--counts", "-c", required=True,
        help="featureCounts output TSV with 'Geneid' and 'Length' columns plus counts"
    )
    parser.add_argument(
        "--bam", "-b", required=True,
        help="Mapped, coordinate-sorted & indexed BAM file"
    )
    parser.add_argument(
        "--out", "-o", required=True,
        help="Output TSV with gene_id, length_bp, raw_counts, RPM, RPKM"
    )
    args = parser.parse_args()

    sys.stderr.write("Counting total mapped reads...\n")
    total = count_mapped_reads(args.bam)
    sys.stderr.write(f"Total mapped reads: {total}\n")

    sys.stderr.write("Parsing featureCounts output...\n")
    lengths, counts = load_featurecounts(args.counts)

    sys.stderr.write("Computing RPM and RPKM...\n")
    results = compute_normalization(lengths, counts, total)

    sys.stderr.write("Writing normalised output...\n")
    write_output(results, args.out)

    sys.stderr.write(f"Done. Results in {args.out}\n")

if __name__ == "__main__":
    main()
