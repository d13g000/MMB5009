"""
clean_msdata.py

This script reads an Excel file containing mass spectrometry (MS) data and
applies the following filters:
  1) The 'Modifications' column must contain at least one occurrence of "{
  n}xMethyl", "{n}xDimethyl" or "{n}xTrimethyl" (for any integer n ≥ 1,
  case-insensitive).
  2) The '# Missed Cleavages' column must be > 0.
  3) The peptide must NOT have any modification on its first or last residue.
  4) The 'Abundance Ratio: (Cancer) / (Healthy)' column must be present,
  parseable as a number, and >= 1.2.
     • If >= 1.5 → the row receives a "High" flag in a new column called
     "Abundance Ratio Flag".
     • If between 1.2 and 1.5 → the row receives a "Medium" flag.
     • Any row with a missing or non-numeric abundance ratio, or abundance ≤
     1.2, is dropped.

Finally, the cleaned DataFrame (only rows meeting all criteria, plus the new
"Abundance Ratio Flag" column) is written to the specified output CSV file.

Usage:
    python3 clean_msdata.py input.xlsx output.csv
    (where input.py is the name of the Excel sheet that is going to be
    filtered and output.csv is the desired name of the CSV file generated.)
"""

import sys
import argparse
import re
import pandas as pd

def has_target_methyl(mod_str: str) -> bool:
    """
    Returns True if the 'Modifications' string contains at least one
    occurrence of:
        {one or more digits}xMethyl
        {one or more digits}xDimethyl
        {one or more digits}xTrimethyl
    (case-insensitive). Otherwise, returns False.

    Any empty or NaN mod_str returns False.
    """
    if pd.isna(mod_str):
        return False

    # Regex: look for “<digits>xMethyl” or “<digits>xDimethyl” or
    # “<digits>xTrimethyl” at a word boundary
    pattern = re.compile(r'\b\d+x(?:Methyl|Dimethyl|Trimethyl)\b',
                         re.IGNORECASE)
    return bool(pattern.search(mod_str))


def has_terminal_mod(annotated_seq: str, mod_str: str) -> bool:
    """
    Returns True if the peptide (extracted from annotated_seq) has any
    modification whose reported position is either 1 or equal to peptide_length.

    If any position == 1 or position == len(peptide), there IS a terminal mod.

    If annotated_seq or mod_str is NaN/None, return False (i.e. “no terminal
    mod detected”).
    """
    if pd.isna(annotated_seq) or pd.isna(mod_str):
        return False

    parts = annotated_seq.split(".")
    if len(parts) < 2:
        # Unexpected format (not “X.Y.Z”), assume no terminal modification
        return False

    peptide = parts[1]
    peptide_length = len(peptide)

    # Regex to capture “[<AA><position>…]”, e.g. “[K1(100)]” or “[M13(99.4)]”
    position_pattern = re.compile(r'\[([A-Za-z])(\d+).*?]', re.IGNORECASE)

    for match in position_pattern.finditer(mod_str):
        pos = int(match.group(2))
        if pos == 1 or pos == peptide_length:
            return True

    return False


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter an MS dataset so that only rows containing at least one of "
            "'<n>xMethyl', '<n>xDimethyl', or '<n>xTrimethyl' in Modifications, "
            "with missed cleavages > 0, no terminal modification, and an "
            "abundance ratio >= 1.2. Rows with abundance >= 1.5 are flagged "
            "'High', 1.2–1.5 flagged 'Medium'."
        )
    )
    parser.add_argument(
        "input_file",
        help="Path to the input Excel file (e.g. data.xlsx)."
    )
    parser.add_argument(
        "output_file",
        help="Path where the filtered CSV file will be saved (e.g. "
             "cleaned_data.csv)."
    )
    args = parser.parse_args()

    # 1) Read the Excel file into a DataFrame
    try:
        df = pd.read_excel(args.input_file, engine="openpyxl")
    except Exception as e:
        print(f"Error reading '{args.input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # 2) Verify required columns exist
    required_cols = {
        "Annotated Sequence",
        "Modifications",
        "# Missed Cleavages",
        "Abundance Ratio: (Cancer) / (Healthy)"
    }
    missing = required_cols - set(df.columns)
    if missing:
        print(f"Input file is missing required columns: "
              f"{', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    # 3) Filter 1 – Must have at least one "<n>xMethyl/Dimethyl/Trimethyl"
    methyl_mask = df["Modifications"].apply(has_target_methyl)
    df_methyl = df[methyl_mask].copy()

    # 4) Filter 2 – "# Missed Cleavages" > 0
    #    Coerce to numeric (invalid parsing → NaN), then keep only > 0
    df_methyl["# Missed Cleavages"] = pd.to_numeric(
        df_methyl["# Missed Cleavages"], errors="coerce"
    )
    missed_mask = df_methyl["# Missed Cleavages"] > 0
    df_cleave = df_methyl[missed_mask].copy()

    # 5) Filter 3 – No terminal modification (pos != 1 and pos !=
    # peptide_length)
    terminal_mask = df_cleave.apply(
        lambda row: not has_terminal_mod(row["Annotated Sequence"],
                                         row["Modifications"]), axis=1
    )
    df_noterm = df_cleave[terminal_mask].copy()

    # 6) Filter 4 – Abundance Ratio > 1.2, drop missing or non-numeric
    #    a) Coerce "Abundance Ratio: (Cancer) / (Healthy)" to numeric
    ratio_col = "Abundance Ratio: (Cancer) / (Healthy)"
    df_noterm["Abundance Ratio"] = pd.to_numeric(
        df_noterm[ratio_col], errors="coerce"
    )

    #    b) Drop any row where Abundance Ratio is NaN (i.e. missing or not
    #    parseable)
    df_ratio = df_noterm.dropna(subset=["Abundance Ratio"]).copy()

    #    c) Keep only ratio >= 1.2
    ratio_mask = df_ratio["Abundance Ratio"] >= 1.2
    df_above12 = df_ratio[ratio_mask].copy()

    # 7) Add a new column "Abundance Ratio Flag": “High” if >= 1.5,
    # else “Medium” (since we know it’s >= 1.2)
    def flag_ratio(x: float) -> str:
        return "High" if x > 1.5 else "Medium"

    df_above12["Abundance Ratio Flag"] = df_above12["Abundance Ratio"].apply(
        flag_ratio)

    # 8) Final DataFrame is df_above12. Write it to CSV.
    try:
        df_above12.to_csv(args.output_file, index=False)
        print(
            f"Filtered data written to '{args.output_file}' "
            f"({len(df_above12)} rows; "
            f"{(df_above12['Abundance Ratio Flag'] == 'High').sum()} High, "
            f"{(df_above12['Abundance Ratio Flag'] == 'Medium').sum()} Medium)."
        )
    except Exception as e:
        print(f"Error writing '{args.output_file}': {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
