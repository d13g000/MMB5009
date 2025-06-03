#!/usr/bin/env python3
"""
annotate_to_fasta.py

Reads a CSV with at least these key columns:
  • "Annotated Sequence"               (e.g. "[K].MDSTEPPYSQKR.[Y]")
  • "Modifications"                    (e.g. "1xDimethyl [K11(99.4)]")
  • "Master Protein Accessions"        (e.g. "P68104")
  • "Positions in Master Proteins"     (e.g. "P68104 [155-166]")
  • "Modifications in Master Proteins" (e.g. "P68104 1xDimethyl [K165(99.4)]")
  • "Abundance Ratio: (Cancer) / (Healthy)" (e.g. 1.45)
  • "Abundance Ratio Flag"             (either "High" or "Medium")

For each row, this script:
  1) Fetches the full UniProt sequence for the given accession.
  2) Parses out the peptide’s start/end (e.g. 155–166).
  3) Parses out every methyl‐type modification (position, type, count).
  4) Builds an output sequence in “pseudo‐color” form:
       – Residues outside the peptide → lowercase
       – Residues inside the peptide → UPPERCASE
       – If a residue inside the peptide is modified, appends "<Me1:count>",
         "<Me2:count>", or "<Me3:count>" immediately after that letter.
       – Ignores modifications that fall outside the peptide window.
  5) Writes one FASTA record per input row, with a header that includes:
       – UniProt accession
       – peptide range
       – number of “[inside]” modifications
       – the numeric abundance ratio
       – the abundance ratio flag (High or Medium)
       – the original UniProt FASTA header

Usage:
    python annotate_to_fasta.py filtered_input.csv output.fasta
"""

import sys
import argparse
import re
import time

import pandas as pd
import requests

# ----------------------------
# Fetch UniProt sequence
# ----------------------------
def fetch_uniprot_sequence(accession: str) -> (str, str):
    """
    Given a UniProt accession (e.g. "P68104"), retrieves the full sequence
    via the UniProt REST API. Returns (header, sequence), where:
      - header is the raw FASTA header line (without the leading '>')
      - sequence is the full amino‐acid string (no line breaks).
    Raises an exception if something goes wrong (e.g. invalid accession or
    network error).
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    response = requests.get(url, timeout=10)
    if response.status_code != 200:
        raise ValueError(f"Could not fetch UniProt entry '{accession}'. "
                         f"HTTP status {response.status_code}.")

    lines = response.text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Unexpected FASTA format for accession "
                         f"{accession!r}.")

    header = lines[0][1:]  # drop the leading '>'
    seq = "".join(lines[1:])  # join the remaining lines
    return header, seq


# ----------------------------
# Parse “Positions in Master Proteins”
# ----------------------------
def parse_peptide_range(pos_string: str, expected_accession: str) -> (int, int):
    """
    The “Positions in Master Proteins” column typically looks like:
      "P68104 [155-166]" or "P11021 [582-596]"

    Verifies that the prefix matches expected_accession, then extracts
    the two numbers inside the brackets, returning (start, end).
    """
    tokens = pos_string.strip().split()
    if len(tokens) != 2:
        raise ValueError(f"Cannot parse peptide range from {pos_string!r}.")
    acc_in_col, range_part = tokens
    if acc_in_col != expected_accession:
        raise ValueError(f"Accession mismatch: expected {expected_accession!r},"
                         f" but got {acc_in_col!r} in '{pos_string}'.")

    m = re.match(r"\[(\d+)-(\d+)\]", range_part)
    if not m:
        raise ValueError(f"Cannot parse bracketed range from {range_part!r}.")
    return int(m.group(1)), int(m.group(2))


# ----------------------------
# Parse “Modifications in Master Proteins”
# ----------------------------
def parse_master_mods(mods_string: str, expected_accession: str) -> dict:
    """
    Example mods_string formats:
      "P68104 1xDimethyl [K165(99.4)]"
      "P14678 2xDimethyl [R108(100); R112(100)]"
      "P11021 1xTrimethyl [K585(100)]"

    Splits into:
      - Prefix (e.g. "P68104") → must equal expected_accession
      - Then parses all occurrences of "<n>x(Methyl|Dimethyl|Trimethyl) [
      <AA><pos>(<score>)]"

    Returns a dict: { position (int) : (tag, count) }, where tag is
    "Me1"/"Me2"/"Me3" corresponding to Methyl/Dimethyl/Trimethyl.
    """
    tokens = mods_string.strip().split(maxsplit=1)
    if len(tokens) != 2:
        raise ValueError(f"Cannot parse Master Mods from {mods_string!r}.")
    acc_in_col, rest = tokens
    if acc_in_col != expected_accession:
        raise ValueError(f"Accession mismatch in Mods column: expected "
                         f"{expected_accession!r}, but got {acc_in_col!r}.")

    # Regex: capture each "(\d+)x(Methyl|Dimethyl|Trimethyl)" followed by "[
    # <AA><pos>(...)]"
    pattern = re.compile(
        r"(\d+)x(Methyl|Dimethyl|Trimethyl)\s*\[([A-Z])(\d+)\([^\]]+\)\]",
        re.IGNORECASE
    )

    result = {}
    for match in pattern.finditer(rest):
        count_str, modname, aa_letter, pos_str = match.groups()
        pos = int(pos_str)
        count = int(count_str)

        m_lower = modname.lower()
        if m_lower == "methyl":
            tag = "Me1"
        elif m_lower == "dimethyl":
            tag = "Me2"
        elif m_lower == "trimethyl":
            tag = "Me3"
        else:
            # This should never happen if the regex is correct
            raise ValueError(f"Unexpected modification type {modname!r} in "
                             f"{mods_string!r}.")

        result[pos] = (tag, count)

    return result


# ----------------------------
# Build “pseudo‐colored” sequence
# ----------------------------
def build_annotated_sequence(
    full_seq: str,
    pep_start: int,
    pep_end: int,
    mods_dict: dict
) -> str:
    """
    full_seq: entire protein sequence (1-based indexing).
    pep_start, pep_end:  1-based start/end of peptide (e.g. 155, 166).
    mods_dict: { position (int) : (tag, count) } e.g. {165: ("Me2",1),
    108:("Me2",2) }.

    Output rules:
      - Residues outside pep_start..pep_end → lowercase
      - Residues inside pep_start..pep_end → UPPERCASE
      - If i in mods_dict AND (pep_start ≤ i ≤ pep_end), append "<{tag}:{
      count}>"
        immediately after the uppercase letter. Mods outside the peptide
        window are ignored.
    """
    annotated = []
    for i, aa in enumerate(full_seq, start=1):
        in_pep = (pep_start <= i <= pep_end)
        base = aa.upper() if in_pep else aa.lower()

        if in_pep and (i in mods_dict):
            tag, count = mods_dict[i]
            annotated.append(f"{base}<{tag}:{count}>")
        else:
            annotated.append(base)

    return "".join(annotated)


# ----------------------------
# Main script
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description=(
            "From a filtered CSV, fetch each Master‐Protein sequence from "
            "UniProt, map the peptide (UPPERCASE) and internal modifications "
            "(inline tags), and write a single FASTA file. Headers include "
            "the abundance ratio and flag."
        )
    )
    parser.add_argument(
        "input_csv",
        help="Path to the filtered CSV (from previous steps)."
    )
    parser.add_argument(
        "output_fasta",
        help="Path where the annotated FASTA will be saved."
    )
    args = parser.parse_args()

    # 1) Read the filtered CSV
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error reading '{args.input_csv}': {e}", file=sys.stderr)
        sys.exit(1)

    # Verify required columns
    required_cols = {
        "Annotated Sequence",
        "Modifications",
        "Master Protein Accessions",
        "Positions in Master Proteins",
        "Modifications in Master Proteins",
        "Abundance Ratio: (Cancer) / (Healthy)",
        "Abundance Ratio Flag"
    }
    missing = required_cols - set(df.columns)
    if missing:
        print(f"Input CSV is missing required columns: {', '.join(missing)}",
              file=sys.stderr)
        sys.exit(1)

    fasta_records = []

    for idx, row in df.iterrows():
        accession = str(row["Master Protein Accessions"]).strip()
        pos_string = str(row["Positions in Master Proteins"]).strip()
        mods_string = str(row["Modifications in Master Proteins"]).strip()

        # Also retrieve abundance ratio and flag
        # If the column is blank or non‐numeric, we force a float NaN,
        # but skip if NaN:
        try:
            ratio = float(row["Abundance Ratio: (Cancer) / (Healthy)"])
        except Exception:
            # If the ratio is missing or non‐numeric, skip this row altogether
            print(f"[Row {idx}] Missing/non‐numeric abundance ratio → "
                  f"skipping.", file=sys.stderr)
            continue

        flag = str(row["Abundance Ratio Flag"]).strip()
        if flag not in {"High", "Medium"}:
            # If flag is something else or missing, we’ll still include it,
            # but note it
            flag = flag or "Unknown"

        # 2) Fetch UniProt sequence (header + sequence)
        try:
            header_line, full_seq = fetch_uniprot_sequence(accession)
        except Exception as e:
            print(f"[Row {idx}] Error fetching UniProt '{accession}': {e}",
                  file=sys.stderr)
            continue

        # 3) Parse peptide range
        try:
            pep_start, pep_end = parse_peptide_range(pos_string, accession)
        except Exception as e:
            print(f"[Row {idx}] Error parsing peptide position '"
                  f"{pos_string}': {e}", file=sys.stderr)
            continue

        # 4) Parse master‐protein modifications
        try:
            mods_dict = parse_master_mods(mods_string, accession)
        except Exception as e:
            print(f"[Row {idx}] Error parsing master modifications "
                  f"'{mods_string}': {e}", file=sys.stderr)
            mods_dict = {}

        # 5) Build the “pseudo‐colored” sequence
        annotated_seq = build_annotated_sequence(full_seq, pep_start,
                                                 pep_end, mods_dict)

        # 6) Build FASTA header. Include:
        #    >accession | peptide=pep_start-pep_end | mods_inside=count_inside |
        #    ratio=<ratio> | flag=<High/Medium> | original_uniprot_header
        count_inside = sum(
            count for pos, (tag, count) in mods_dict.items()
            if pep_start <= pos <= pep_end
        )

        header = (
            f"{accession} "
            f"| peptide={pep_start}-{pep_end} "
            f"| mods_inside={count_inside} "
            f"| ratio={ratio:.3f} "
            f"| flag={flag} "
            f"| {header_line}"
        )

        fasta_records.append((header, annotated_seq))

        # Polite pause so we don’t overwhelm UniProt
        time.sleep(0.5)

    # 7) Write out all FASTA records (with a blank line between entries)
    try:
        with open(args.output_fasta, "w") as fout:
            for header, seq in fasta_records:
                fout.write(f">{header}\n")
                # wrap at 60 characters per line for readability
                for i in range(0, len(seq), 60):
                    fout.write(seq[i : i + 60] + "\n")
                fout.write("\n")  # extra blank line between entries
    except Exception as e:
        print(f"Error writing '{args.output_fasta}': {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote {len(fasta_records)} FASTA entries to '{args.output_fasta}'.")


if __name__ == "__main__":
    main()