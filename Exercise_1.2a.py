import sys  # Import the sys module to handle command-line arguments


def trim_sequence(sequence, mode, n):
    """Trim the sequence based on mode and number of bases."""
    if mode == "start":  # If trimming from the start
        trimmed_seq = sequence[n:]  # Remove the first n bases
    elif mode == "end":  # If trimming from the end
        trimmed_seq = sequence[
                      :-n] if n > 0 else sequence  # Remove the last n bases
    elif mode == "both":  # If trimming from both ends
        trimmed_seq = sequence[
                      n:-n] if n > 0 else sequence  # Remove n bases from both ends
    else:
        print(
            "Invalid mode. Use 'start', 'end', or 'both'.")  # Handle invalid mode input
        sys.exit(1)  # Exit the script with an error

    return trimmed_seq  # Return the trimmed sequence


def process_fastq(file_path, mode, n):
    """Read and process the FASTQ file."""
    try:
        with open(file_path,
                  'r') as file:  # Open the input FASTQ file for reading
            lines = file.readlines()  # Read all lines from the file

        output_lines = []  # List to store the modified FASTQ lines
        for i in range(0, len(lines),
                       4):  # Loop through the file in FASTQ format (4 lines per entry)
            header = lines[i].strip()  # First line: Sequence identifier
            sequence = lines[i + 1].strip()  # Second line: DNA sequence
            plus = lines[i + 2].strip()  # Third line: Separator (+)
            quality = lines[i + 3].strip()  # Fourth line: Quality scores

            trimmed_seq = trim_sequence(sequence, mode,
                                        n)  # Trim the sequence based on the given mode
            trimmed_qual = trim_sequence(quality, mode,
                                         n)  # Trim the quality scores similarly

            if not trimmed_seq:  # If trimming removes the entire sequence
                print(
                    f"Warning: Sequence at {header} has been completely removed.")  # Warn the user
                continue  # Skip writing this entry to the output file

            output_lines.extend([header, trimmed_seq, plus,
                                 trimmed_qual])  # Add trimmed entry to output list

        output_file = file_path.replace(".fastq",
                                        "_trimmed.fastq")  # Create output file name
        with open(output_file, 'w') as file:  # Open the output file for writing
            file.write(
                "\n".join(output_lines) + "\n")  # Write the modified FASTQ data

        print(
            f"Trimming complete. Output saved to {output_file}")  # Notify user of success
    except Exception as e:
        print(f"Error: {e}")  # Print any errors that occur during processing


if __name__ == "__main__":  # Ensure script runs only when executed directly
    if len(sys.argv) != 4:  # Check if the correct number of arguments are provided
        print(
            "Usage: python script.py <input_file> <mode> <n>")  # Print usage instructions
        print("<mode>: 'start', 'end', or 'both'")  # Explain mode options
        sys.exit(1)  # Exit script with an error

    input_file = sys.argv[
        1]  # Get the input FASTQ file name from command-line argument
    mode = sys.argv[2].lower()  # Get the trimming mode and convert to lowercase
    n = int(
        sys.argv[3])  # Get the number of bases to trim, converted to an integer

    process_fastq(input_file, mode,
                  n)  # Call function to process the FASTQ file
