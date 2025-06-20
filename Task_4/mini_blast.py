"""
mini_blast.py

A simplified, BLAST-like search for protein sequences.
Follows Pertsemlidis et al. (2001) “Having a BLAST with bioinformatics (and
avoiding BLASTphemy)”. Steps implemented:
  1) Read input FASTA files
  2) Define a substitution scoring matrix (BLOSUM62)
  3) Seed: scan for all exact k-mer matches
  4) Extend each seed in both directions, accumulate score
  5) Record the best High‐scoring Segment Pair (HSP) per target
  6) Rank targets by HSP score
  7) Output ranking and print query header (to see its annotation/function)

Simplifying assumptions:
  - Configurable k-mer/word size (-k / --wordsize, default=3) required for exact
  matching to seed.
  - Choice of ungapped (default) or gapped alignment (--ungapped / --gapped)
  - Uses full BLOSUM62 for scoring matches/mismatches.
  - Threshold for reporting any HSP is > 0.
"""

import argparse
from collections import defaultdict
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices

# -------- Step 1: Read sequences ----------
def read_fasta(filename):
    """
    Reads one-sequence FASTA and returns (header, seq) or list of
    SeqRecord.
    """
    return list(SeqIO.parse(filename, "fasta"))

# -------- Step 2: Load scoring matrix ----------
blosum62 = substitution_matrices.load("BLOSUM62")

def score_pair(a, b):
    """
    Returns BLOSUM62 score for a pair of residues (or -4 if undefined).
    """
    key = (a, b)
    if key not in blosum62:
        key = (b, a)
    return blosum62.get(key, -4)

# -------- Step 3 & 4: Seed and extend ----------
def find_kmer_positions(seq, k):
    """
    Returns a dictionary mapping each k-mer to list of start positions in the
    sequence.
    """
    d = defaultdict(list)
    for i in range(len(seq) - k + 1):
        word = str(seq[i:i+k])
        d[word].append(i)
    return d

def extend_hsp_ungapped(query, target, qpos, tpos, k):
    """
    Ungapped extension of a seed at (qpos, tpos) of length k.

    X-drop style: stops extending when cumulative score < 0 and returns the
    best score found.
    """
    # initial k-mer score
    cur = sum(score_pair(query[qpos + i], target[tpos + i]) for i in range(k))
    best = cur

    # extend right
    i = 1
    while qpos + k - 1 + i < len(query) and tpos + k - 1 + i < len(target):
        cur += score_pair(query[qpos + k - 1 + i], target[tpos + k - 1 + i])
        if cur < 0:
            break
        best = max(best, cur)
        i += 1

    # extend left (starting from best right-extended score)
    cur = best
    i = 1
    while qpos - i >= 0 and tpos - i >= 0:
        cur += score_pair(query[qpos - i], target[tpos - i])
        if cur < 0:
            break
        best = max(best, cur)
        i += 1

    return best

def ungapped_blast(query_seq, target_seq, k):
    """
    Computes the best ungapped HSP score between query and target.
    """
    best_score = 0
    tgt_index = find_kmer_positions(target_seq, k)
    for i in range(len(query_seq) - k + 1):
        word = str(query_seq[i : i + k])
        for tpos in tgt_index.get(word, []):
            sc = extend_hsp_ungapped(query_seq, target_seq, i, tpos, k)
            if sc > best_score:
                best_score = sc
    return best_score

def gapped_blast(query_seq, target_seq):
    """
    Computes the best gapped, local alignment score using BLOSUM62 and affine
    gaps where:
    - gap open = -10
    - gap extend = -1
    """
    alignments = pairwise2.align.localds(
        query_seq,
        target_seq,
        blosum62,
        -10,   # gap open penalty
        -1,    # gap extend penalty
        one_alignment_only=True,
        score_only=True
    )
    # if no alignment found, score_only returns 0
    return alignments if alignments is not None else 0

def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("query", help="Query FASTA (one sequence)")
    p.add_argument("target", help="Target FASTA (multiple "
                                  "sequences)")
    p.add_argument(
        "-k", "--wordsize",
        type=int,
        default=3,
        help="Word size for seeding (only used for ungapped)"
    )
    grp = p.add_mutually_exclusive_group()
    grp.add_argument(
        "--gapped",
        dest="gapped",
        action="store_true",
        help="Use gapped local alignment"
    )
    grp.add_argument(
        "--ungapped",
        dest="gapped",
        action="store_false",
        help="Use ungapped HSP extension"
    )
    p.set_defaults(gapped=False)  # default = ungapped

    args = p.parse_args()
    # Read sequences
    qrec = read_fasta(args.query)[0]
    targets = read_fasta(args.target)

    # Choose alignment function
    if args.gapped:
        align_fn = gapped_blast
        mode = "GAPPED"
    else:
        align_fn = lambda q, t: ungapped_blast(q, t, args.wordsize)
        mode = f"UNGAPPED (k={args.wordsize})"

    # Score all targets
    results = []
    for trec in targets:
        sc = align_fn(str(qrec.seq), str(trec.seq))
        results.append((trec.id, trec.description, sc))

    # Sort descending by score
    results.sort(key=lambda x: x[2], reverse=True)

    # Output
    print(f"Query: {qrec.id}")
    print(f"Description: {qrec.description}")
    print(f"Mode: {mode}\n")

    print("Rank\tTarget_ID\tScore\tDescription")
    for i, (tid, desc, sc) in enumerate(results, 1):
        print(f"{i}\t{tid}\t{sc:.1f}\t{desc}")

if __name__ == "__main__":
    main()