"""
bm_search.py

Usage:
    python3 bm_search.py <reference_fasta> <pattern> [--mismatches N] [
    --benchmark]

Description:
    - Loads a reference genome FASTA file into memory.
    - Normalizes text and pattern to uppercase.
    - Implements Boyer–Moore exact matching with Bad-Character and
    Good-Suffix heuristics, always choosing the larger of the two shifts on
    mismatch, and counts how many alignments
      it actually inspects (to compare against naive).
    - Implements a naive sliding-window matcher (counts its alignments too).
    - Wraps Boyer–Moore in a Pigeon-Hole strategy to allow up to N mismatches.
    - Benchmarks naive, Boyer–Moore, and Pigeon-Hole (when mismatches > 0).
    - Reports number of hits for each method, how many alignments BM
    skipped, and makes runtimes >60s human-readable.
"""

import sys
import time
import argparse

def load_fasta(filepath):
    """Read FASTA, concatenate sequence lines, uppercase."""
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith('>'):
                seq.append(line.strip().upper())
    return ''.join(seq)

############################
# Boyer-Moore preprocessing
############################

def build_bad_char_table(pat):
    table = {chr(i): -1 for i in range(256)}
    for i, ch in enumerate(pat):
        table[ch] = i
    return table

def build_good_suffix_tables(pat):
    m = len(pat)
    suff = [0] * m
    suff[m-1] = m
    g = f = m-1
    for i in range(m-2, -1, -1):
        if i > g and suff[i + m-1-f] < i-g:
            suff[i] = suff[i + m-1-f]
        else:
            if i < g:
                g = i
            f = i
            while g >= 0 and pat[g] == pat[g + m-1-f]:
                g -= 1
            suff[i] = f-g
    good_suffix = [m] * m
    border_shift = [0] * (m+1)
    j = 0
    for i in range(m-1, -1, -1):
        if suff[i] == i+1:
            while j < m-1-i:
                if border_shift[j] == 0:
                    border_shift[j] = m-1-i
                j += 1
    for i in range(m-1):
        good_suffix[m-1-suff[i]] = m-1-i
    for i in range(m+1):
        if border_shift[i] == 0:
            border_shift[i] = m
    return good_suffix, border_shift

def boyer_moore_search(text, pat):
    """
    Returns (matches, alignments_inspected):
      - matches: list of match positions
      - alignments_inspected: how many shifts s were actually tested
    """
    n, m = len(text), len(pat)
    if m == 0:
        return list(range(n+1)), n+1

    bad_char = build_bad_char_table(pat)
    good_suffix, border_shift = build_good_suffix_tables(pat)

    matches = []
    s = 0
    alignments = 0
    while s <= n - m:
        alignments += 1
        j = m - 1
        while j >= 0 and pat[j] == text[s + j]:
            j -= 1
        if j < 0:
            matches.append(s)
            s += border_shift[0]
        else:
            bc_shift = j - bad_char.get(text[s+j], -1)
            gs_shift = good_suffix[j]
            s += max(bc_shift, gs_shift)
    return matches, alignments

##########################
# Naive exact matching
##########################

def naive_search(text, pat):
    """
    Returns (matches, alignments_inspected):
      - matches: list of match positions
      - alignments_inspected: always n-m+1 (or n+1 if m==0)
    """
    n, m = len(text), len(pat)
    if m == 0:
        return list(range(n+1)), n+1
    hits = []
    alignments = 0
    for i in range(n - m + 1):
        alignments += 1
        if text[i:i+m] == pat:
            hits.append(i)
    return hits, alignments

#######################################
# Pigeon-Hole Principle for n mismatches
#######################################

def pigeonhole_search(text, pat, max_mm):
    n, m = len(text), len(pat)
    if m == 0:
        return [(i,0) for i in range(n+1)]
    k = max_mm
    pieces = k + 1
    q, r = divmod(m, pieces)
    splits = []
    start = 0
    for i in range(pieces):
        length = q + (1 if i < r else 0)
        splits.append((start, start + length))
        start += length

    candidates = set()
    for a, b in splits:
        sub = pat[a:b]
        if not sub:
            continue
        hits, _ = boyer_moore_search(text, sub)
        for h in hits:
            align = h - a
            if 0 <= align <= n - m:
                candidates.add(align)

    results = []
    for pos in sorted(candidates):
        window = text[pos:pos+m]
        mismatches = sum(1 for x,y in zip(window, pat) if x != y)
        if mismatches <= k:
            results.append((pos, mismatches))
    return results

#############################
# Timing / Performance Utils
#############################

def format_time(t):
    """If t>60s, render as Hh Mm Ss; else as 'X.XXXX s'."""
    if t > 60:
        hrs = int(t // 3600)
        rem = t % 3600
        mins = int(rem // 60)
        secs = rem % 60
        parts = []
        if hrs:
            parts.append(f"{hrs}h")
        if mins or hrs:
            parts.append(f"{mins}m")
        parts.append(f"{secs:.3f}s")
        return " ".join(parts)
    else:
        return f"{t:.4f} s"

def timeit(fn, *args):
    t0 = time.perf_counter()
    res = fn(*args)
    t1 = time.perf_counter()
    return res, t1 - t0

###########################
# Command-Line Interface
###########################

def main():
    parser = argparse.ArgumentParser(
        description="Boyer–Moore vs Naive Search with Pigeon-Hole  "
                    "approximate matching"
    )
    parser.add_argument('reference_fasta', help='Path to reference '
                                                'FASTA')
    parser.add_argument('pattern',       help='Pattern to search '
                                              'for')
    parser.add_argument('--mismatches', '-m', type=int, default=0,
                        help='Allow up to N mismatches')
    parser.add_argument('--benchmark', '-b', action='store_true',
                        help='Benchmark naive, BM, and Pigeon-Hole (if m>0)')
    args = parser.parse_args()

    text = load_fasta(args.reference_fasta)
    pat  = args.pattern.strip().upper()
    n, m = len(text), len(pat)
    print(f"Reference length: {n} bp", file=sys.stderr)
    print(f"Pattern length:   {m} bp", file=sys.stderr)

    # Exact matching
    if args.mismatches == 0:
        print("\n=== Exact Matching ===")
        naive_res, t_naive = timeit(naive_search, text, pat)
        bm_res,    t_bm    = timeit(boyer_moore_search, text, pat)

        naive_hits, naive_align = naive_res
        bm_hits,    bm_align    = bm_res

        print(f"Naive:        {len(naive_hits)} hits, inspected "
              f"{naive_align} alignments, time = {format_time(t_naive)}")
        print(f"Boyer–Moore:  {len(bm_hits)} hits, inspected {bm_align} "
              f"alignments, time = {format_time(t_bm)}")
        skipped = naive_align - bm_align
        print(f"BM skipped {skipped} alignments compared to naive")

        if args.benchmark:
            if set(naive_hits) != set(bm_hits):
                print("WARNING: result sets differ!", file=sys.stderr)
            else:
                print("Results match exactly.", file=sys.stderr)

        print(f"\nMatches at positions ({len(bm_hits)} total):")
        print("\n".join(map(str, bm_hits)))

    # Approximate matching
    else:
        k = args.mismatches
        print(f"\n=== Approximate Matching (≤ {k} mismatches) ===")
        # always run timings
        naive_res, t_naive = timeit(naive_search, text, pat)
        bm_res,    t_bm    = timeit(boyer_moore_search, text, pat)
        ph_res,    t_ph    = timeit(pigeonhole_search, text, pat, k)

        naive_hits, naive_align = naive_res
        bm_hits,    bm_align    = bm_res
        ph_hits                = ph_res

        print(f"Naive (exact):        {len(naive_hits)} hits, inspected "
              f"{naive_align} alignments, time = {format_time(t_naive)}")
        print(f"Boyer–Moore (exact):  {len(bm_hits)} hits, inspected "
              f"{bm_align} alignments, time = {format_time(t_bm)}")
        skipped = naive_align - bm_align
        print(f"\nBM skipped {skipped} alignments compared to naive")
        print(f"\nPigeon-Hole (≤{k} mm): {len(ph_hits)} hits, time = "
              f"{format_time(t_ph)}")

        print(f"\nFound {len(ph_hits)} hits (pos → mismatches):")
        for pos, mm in ph_hits:
            print(f"{pos}\t{mm}")

if __name__ == '__main__':
    main()
