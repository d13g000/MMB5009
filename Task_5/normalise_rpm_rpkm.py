#!/usr/bin/env python3
import argparse
import sys
import pysam

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute RPM and RPKM from featureCounts output."
    )
    p.add_argument(
        "--counts", "-c", required=True,
        help="featureCounts output (tab) with columns: Geneid,Chr,Start,End,"
             "Strand,Length,<sample>"
    )
    p.add_argument(
        "--bam", "-b", required=True,
        help="Mapped BAM file (coordinate-sorted & indexed)"
    )
    p.add_argument(
        "--out", "-o", required=True,
        help="Output TSV with gene_id, length_bp, raw_counts, RPM, RPKM"
    )
    return p.parse_args()

def count_mapped_reads(bam_path):
    """Return total mapped reads (excludes unmapped, secondary, qcfail,
    duplicates)."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    total = 0
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_qcfail or r.is_duplicate:
            continue
        total += 1
    bam.close()
    return total

def load_featurecounts(counts_path):
    """
    Parse featureCounts output.
    Returns:
      gene_lengths: dict[gene_id] = length_bp (int)
      raw_counts:   dict[gene_id] = count (int)
    """
    gene_lengths = {}
    raw_counts   = {}
    with open(counts_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            # header line
            if parts[0] == "Geneid":
                # find indices
                idx_len   = parts.index("Length")
                idx_count = len(parts) - 1
                continue
            gene = parts[0]
            length = int(parts[idx_len])
            count  = int(parts[idx_count])
            gene_lengths[gene] = length
            raw_counts[gene]   = count
    return gene_lengths, raw_counts

def write_norm(gene_lengths, raw_counts, total_mapped, out_path):
    """Compute RPM/RPKM and write to TSV."""
    million = total_mapped / 1e6
    with open(out_path, "w") as out:
        out.write("gene_id\tlength_bp\traw_counts\tRPM\tRPKM\n")
        for gene, length in gene_lengths.items():
            cnt = raw_counts.get(gene, 0)
            rpm  = cnt / million if million > 0 else 0
            rpkm = rpm / (length/1000) if length > 0 else 0
            out.write(f"{gene}\t{length}\t{cnt}\t{rpm:.6f}\t{rpkm:.6f}\n")

def main():
    args = parse_args()
    sys.stderr.write("Counting total mapped reads...\n")
    total = count_mapped_reads(args.bam)
    sys.stderr.write(f"Total mapped reads: {total}\n")
    sys.stderr.write("Parsing featureCounts output...\n")
    lengths, counts = load_featurecounts(args.counts)
    sys.stderr.write("Computing RPM/RPKM and writing output...\n")
    write_norm(lengths, counts, total, args.out)
    sys.stderr.write(f"Done. Results in {args.out}\n")

if __name__ == "__main__":
    main()
