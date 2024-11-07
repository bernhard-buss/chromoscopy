def split_gtf_by_chromosome(gtf_filename):
    """
    Splits a GTF file by chromosome and writes each chromosome's features into separate files.

    Parameters:
    - gtf_filename (str): Path to the input GTF file.

    Output:
    - For each chromosome, a separate GTF file containing all features for that chromosome.
    """
    import os

    # Dictionary to hold open file handles for each chromosome
    chrom_files = {}
    # List to hold header lines (lines starting with '#')
    header_lines = []

    try:
        with open(gtf_filename, 'r') as infile:
            for line in infile:
                # Check if the line is a header (starts with '#')
                if line.startswith('#'):
                    header_lines.append(line)
                    continue  # Skip to the next line

                # Split the line into fields
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    # Not a valid GTF line, skip it
                    continue

                # Extract the chromosome name from the first column
                chrom = fields[0]

                # Ensure we have a file open for this chromosome
                if chrom not in chrom_files:
                    # Generate a unique filename for the chromosome
                    base_name = os.path.splitext(gtf_filename)[0]
                    output_filename = f"{base_name}_{chrom}.gtf"
                    # Open a new file for writing
                    chrom_files[chrom] = open(output_filename, 'w')
                    # Write header lines to the new chromosome file
                    chrom_files[chrom].writelines(header_lines)

                # Write the current line to the chromosome's file
                chrom_files[chrom].write(line)

    finally:
        # Close all open file handles
        for f in chrom_files.values():
            f.close()
