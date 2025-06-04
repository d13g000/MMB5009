#!/usr/bin/env python3
"""
generate_annotated_csv.py

Reads a filtered MS‐data CSV and produces a new CSV with these columns:
  - Peptide: a multiline field with
      Control :
      <control sequence>

      Cancer :
      <cancer sequence>
  - Modification: copied directly from the input "Modifications" column
  - Positions in Master Proteins: copied directly from input
  - Miscleavage: copied from “# Missed Cleavages”
  - Sample: a multiline semicolon-separated list of sample IDs;
      Control:
      [S1]; [S2]; [S3];

      Cancer:
      [S4]; [S5]; [S6]
      (only including samples if value != "Not Found"
  - Abundance: a multiline field with
      Control:
      <Abundances (Grouped): Healthy>

      Cancer:
      <Abundances (Grouped): Cancer>
  - PSMs: copied from “# PSMs”

Usage:
    python3 generate_annotated_csv.py input.csv output.csv
    (where input.py is the name of the CSV file that is going to be read [
    produced by clean_msdata.py script] and output.csv is the desired
    name of the annotated CSV file generated.)
"""

import sys
import argparse
import re
import pandas as pd

def extract_peptide(annotated_seq: str) -> str:
    """
    Return peptide sequence from "Annotated Sequence" column.
    """
    if pd.isna(annotated_seq):
        return ""
    parts = annotated_seq.split(".")
    if len(parts) < 2:
        return annotated_seq.strip()
    return parts[1].strip()

def parse_methyl_mods(mod_str: str):
    """
    Parse only methyl‐type modifications from the “Modifications” column.
    """
    if pd.isna(mod_str):
        return []

    result = []
    block_pattern = re.compile(
        r"(\d+)x(Methyl|Dimethyl|Trimethyl)\s*\[([^]]+)]",
        re.IGNORECASE
    )
    for block_match in block_pattern.finditer(mod_str):
        count_str, modname, inside = block_match.groups()
        mod_lc = modname.lower()
        if mod_lc == "methyl":
            tag = "Me1"
        elif mod_lc == "dimethyl":
            tag = "Me2"
        else:  # "trimethyl"
            tag = "Me3"
        # Inside may contain multiple positions separated by ";" or commas →
        # find all digits before "("
        pos_pattern = re.compile(r"[A-Z](\d+)\(", re.IGNORECASE)
        for pos_match in pos_pattern.finditer(inside):
            pos = int(pos_match.group(1))
            result.append((tag, pos))
    return result

def build_cancer_peptide(base_pep: str, methyl_mods):
    """
    Insert "<{tag}>" immediately after the amino acid at each 1-based
    position 'pos' returning the modified peptide string.
    """
    if not methyl_mods:
        return base_pep

    chars = list(base_pep)
    # Sort by descending position so insertions don't shift later indices
    for tag, pos in sorted(methyl_mods, key=lambda x: x[1], reverse=True):
        idx = pos - 1
        if 0 <= idx < len(chars):
            chars[idx] = f"{chars[idx]}<{tag}>"
    return "".join(chars)

def extract_methyl_substrings(mod_str: str) -> str:
    """
    Extract and rejoin only methylation‐type substrings.
    """
    if pd.isna(mod_str):
        return ""
    pattern = re.compile(
        r"(\d+x(?:Methyl|Dimethyl|Trimethyl)\s*\[[^]]+])",
        re.IGNORECASE
    )
    matches = pattern.findall(mod_str)
    return "; ".join(matches)

def combine_sample_columns(row):
    """
    Build multiline sample field:
      Control:
      [S1]; [S2]; [S3];

      Cancer:
      [S4]; [S5]; [S6]
    Only include sample the corresponding "Found in Sample:" cell is not "Not
    Found".
    """
    control_ids = []
    cancer_ids = []
    for col in row.index:
        if col.startswith("Found in Sample:"):
            val = row[col]
            if pd.notna(val) and str(val).strip() != "Not Found":
                m = re.search(r"\[S(\d+)]", col)
                if m:
                    idx = int(m.group(1))
                    sample_label = f"[S{idx}]"
                    if 1 <= idx <= 3:
                        control_ids.append(sample_label)
                    else:
                        cancer_ids.append(sample_label)

    control_str = "; ".join(control_ids) + ("; " if control_ids else "")
    cancer_str = "; ".join(cancer_ids) + ("; " if cancer_ids else "")

    return (
        "Control:\n"
        f"{control_str}\n\n"
        "Cancer:\n"
        f"{cancer_str}"
    )

def main():
    parser = argparse.ArgumentParser(
        description="Generate a combined control/cancer peptide CSV from "
                    "filtered MS data."
    )
    parser.add_argument("input_csv", help="Path to the filtered "
                                          "input CSV.")
    parser.add_argument("output_csv", help="Path to the output "
                                           "CSV to create.")
    args = parser.parse_args()

    # Read the input CSV
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error reading '{args.input_csv}': {e}", file=sys.stderr)
        sys.exit(1)

    # Verify required columns exist
    required_cols = {
        "Annotated Sequence",
        "Modifications",
        "# Missed Cleavages",
        "Positions in Master Proteins",
        "Abundances (Grouped): Healthy",
        "Abundances (Grouped): Cancer",
        "# PSMs"
    }
    missing = required_cols - set(df.columns)
    if missing:
        print(f"Input CSV is missing required columns: {', '.join(missing)}",
              file=sys.stderr)
        sys.exit(1)

    output_rows = []
    for _, row in df.iterrows():
        base_pep = extract_peptide(row["Annotated Sequence"])
        methyl_mods = parse_methyl_mods(row["Modifications"])

        control_pep = base_pep
        cancer_pep = build_cancer_peptide(base_pep, methyl_mods)

        # Build the "Peptide" multiline field
        peptide_field = (
            "Control :\n"
            f"{control_pep}\n\n"
            "Cancer :\n"
            f"{cancer_pep}"
        )

        modification_field = extract_methyl_substrings(row["Modifications"])
        positions_in_master = row["Positions in Master Proteins"]
        miscleavage = row["# Missed Cleavages"]

        sample_field = combine_sample_columns(row)

        abundance_healthy = row["Abundances (Grouped): Healthy"]
        abundance_cancer = row["Abundances (Grouped): Cancer"]
        # Build the "Abundance" multiline field
        abundance_field = (
            "Control:\n"
            f"{abundance_healthy}\n\n"
            "Cancer:\n"
            f"{abundance_cancer}"
        )

        psms = row["# PSMs"]

        output_rows.append({
            "Control/Cancer Sample Peptide": peptide_field,
            "Modification (and which amino acid is methylated + Probability)":
                modification_field,
            "Positions in Master Proteins": positions_in_master,
            "Miscleavage": miscleavage,
            "Which Sample Replicate it was identified in": sample_field,
            "Abundance": abundance_field,
            "PSMs": psms
        })

    out_df = pd.DataFrame(output_rows)
    out_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(out_df)} rows to '{args.output_csv}'.")

if __name__ == "__main__":
    main()
