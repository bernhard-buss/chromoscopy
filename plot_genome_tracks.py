# plot_genome_tracks.py

# NOTE: implementation not working!

import tempfile
import subprocess
import os
from collections import defaultdict
import configparser

def plot_genome_tracks(processed_data, output_dir="chromosome_plots", genome="hg38", track_order=None, figsize=(10, 6)):
    """
    Plots chromosomes with annotations using pyGenomeTracks for each chromosome individually.

    Parameters:
    - processed_data (dict): Dictionary containing processed genomic data with the following keys:
        - 'exons': List of dictionaries with keys ['chromosome', 'strand', 'gene_id', 'gene_name', 'transcript_id', 'start', 'end']
        - 'introns': List of dictionaries with keys ['chromosome', 'strand', 'gene_id', 'gene_name', 'transcript_id', 'start', 'end']
        - 'intergenic': List of dictionaries with keys ['chromosome', 'strand', 'start', 'end']
        - 'gene_annotations': List of dictionaries with keys ['chromosome', 'strand', 'gene_id', 'gene_name', 'start', 'end', 'midpoint']
    - output_dir (str): Directory to save the generated plot images.
    - genome (str): Genome assembly identifier (e.g., 'hg38', 'mm10'). Ensure it's supported by pyGenomeTracks.
    - track_order (list of str): Optional. Order of tracks to display. Defaults to ['intergenic', 'genes', 'introns', 'exons'].
    - figsize (tuple): Size of the output image in inches (width, height).

    Returns:
    - None. Saves the plots to the specified output directory.
    """

    # Set default track order if not provided
    if track_order is None:
        track_order = ['intergenic', 'genes', 'introns', 'exons']

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Extract unique chromosomes
    chromosomes = set()
    for feature in processed_data.get('exons', []):
        chromosomes.add(feature['chromosome'])
    for feature in processed_data.get('introns', []):
        chromosomes.add(feature['chromosome'])
    for feature in processed_data.get('intergenic', []):
        chromosomes.add(feature['chromosome'])
    for feature in processed_data.get('gene_annotations', []):
        chromosomes.add(feature['chromosome'])

    chromosomes = sorted(chromosomes)  # You can use natsorted if needed

    # Organize data by chromosome
    data_by_chrom = defaultdict(lambda: defaultdict(list))
    for feature in processed_data.get('exons', []):
        data_by_chrom[feature['chromosome']]['exons'].append(feature)
    for feature in processed_data.get('introns', []):
        data_by_chrom[feature['chromosome']]['introns'].append(feature)
    for feature in processed_data.get('intergenic', []):
        data_by_chrom[feature['chromosome']]['intergenic'].append(feature)
    for feature in processed_data.get('gene_annotations', []):
        data_by_chrom[feature['chromosome']]['gene_annotations'].append(feature)

    # Iterate over each chromosome and plot
    for chrom in chromosomes:
        chrom_data = data_by_chrom[chrom]
        
        # Determine the plotting region based on gene annotations
        gene_annotations = chrom_data.get('gene_annotations', [])
        if not gene_annotations:
            print(f"No gene annotations found for {chrom}. Skipping.")
            continue

        # Ensure that all gene annotations have start < end
        valid_gene_annotations = [gene for gene in gene_annotations if gene['start'] < gene['end']]
        if not valid_gene_annotations:
            print(f"No valid gene annotations (start < end) found for {chrom}. Skipping.")
            continue

        max_end = max(feature['end'] for feature in valid_gene_annotations)
        min_start = min(feature['start'] for feature in valid_gene_annotations)
        region = f"{chrom}:{min_start}-{max_end}"

        # Define output file name
        output_file = os.path.join(output_dir, f"{chrom}_plot.png")

        # Create a temporary directory to store BED files and config
        with tempfile.TemporaryDirectory() as temp_dir:
            # Define BED file paths
            bed_files = {
                'exons': os.path.join(temp_dir, 'exons.bed'),
                'introns': os.path.join(temp_dir, 'introns.bed'),
                'intergenic': os.path.join(temp_dir, 'intergenic.bed'),
                'genes': os.path.join(temp_dir, 'genes.bed')
            }

            # Helper function to write BED files
            def write_bed_file(data_list, bed_path, feature_type):
                """
                Writes a BED file from a list of dictionaries.

                Parameters:
                - data_list (list of dict): List containing genomic features.
                - bed_path (str): Path to the output BED file.
                - feature_type (str): Type of feature ('exon', 'intron', 'intergenic', 'gene').

                Returns:
                - None
                """
                with open(bed_path, 'w') as bed_file:
                    for feature in data_list:
                        chrom_feat = feature['chromosome']
                        if chrom_feat != chrom:
                            continue  # Ensure features are for the current chromosome
                        start = int(feature['start']) - 1  # BED format is 0-based
                        end = int(feature['end'])

                        if feature_type == 'intergenic':
                            name = 'intergenic'
                            strand = '.'  # Intergenic regions typically have no strand
                        elif feature_type == 'gene':
                            name = feature.get('gene_name', feature.get('gene_id', 'gene'))
                            strand = feature.get('strand', '.')
                        else:
                            # For exons and introns
                            name = feature.get('gene_name', feature.get('gene_id', 'feature'))
                            strand = feature.get('strand', '.')

                        score = '.'  # Placeholder as score is not provided
                        bed_line = f"{chrom_feat}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
                        bed_file.write(bed_line)

            # Write each BED file
            write_bed_file(chrom_data.get('intergenic', []), bed_files['intergenic'], 'intergenic')
            write_bed_file(chrom_data.get('gene_annotations', []), bed_files['genes'], 'gene')
            write_bed_file(chrom_data.get('introns', []), bed_files['introns'], 'intron')
            write_bed_file(chrom_data.get('exons', []), bed_files['exons'], 'exon')

            # Create the pyGenomeTracks configuration file with proper sections using configparser
            config = configparser.ConfigParser()
            config['global'] = {
                'genome': genome,
                'outFileName': output_file,
                'outFileFormat': 'png',
                'dpi': '300',
                'tracks': ",".join(track_order)
            }

            # Define track configurations
            track_configs = {
                'intergenic': {
                    'file': bed_files['intergenic'],
                    'color': 'lightgrey',
                    'name': 'Intergenic',
                    'type': 'bed',
                    'height': '10'
                },
                'genes': {
                    'file': bed_files['genes'],
                    'color': 'seagreen',
                    'name': 'Genes',
                    'type': 'bed',
                    'itemRgb': 'On',
                    'height': '15'
                },
                'introns': {
                    'file': bed_files['introns'],
                    'color': 'lightblue',
                    'name': 'Introns',
                    'type': 'bed',
                    'height': '10'
                },
                'exons': {
                    'file': bed_files['exons'],
                    'color': 'blue',
                    'name': 'Exons',
                    'type': 'bed',
                    'height': '10'
                }
            }

            for track in track_order:
                if track not in track_configs:
                    print(f"Warning: Track '{track}' not found in track_configs.")
                    continue
                config[track] = {
                    'file': track_configs[track]['file'],
                    'color': track_configs[track]['color'],
                    'name': track_configs[track]['name'],
                    'type': track_configs[track]['type'],
                    'height': track_configs[track]['height']
                }
                if 'itemRgb' in track_configs[track]:
                    config[track]['itemRgb'] = track_configs[track]['itemRgb']

            # Write the configuration file
            config_path = os.path.join(temp_dir, 'pygenometracks.ini')
            with open(config_path, 'w') as config_file:
                config.write(config_file)

            # Prepare the pyGenomeTracks command
            cmd = [
                'pyGenomeTracks',
                '--tracks', config_path,
                '--region', region,
                '-o', output_file
            ]

            # Execute the command
            try:
                result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"Chromosome plot for {chrom} saved to {output_file}")
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while generating the chromosome plot for {chrom} with pyGenomeTracks.")
                print("Standard Output:", e.stdout.decode())
                print("Standard Error:", e.stderr.decode())
                raise
