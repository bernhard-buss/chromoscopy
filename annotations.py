import os
import pickle
import pandas as pd
from natsort import natsorted
import re
import load_gtf

# GTF Base filename pattern
gtf_base_filename = './gtf/Homo_sapiens.GRCh38.113'

def get_annotations_filename(chrom):
    return f"{gtf_base_filename}_{chrom}.gtf"

def preprocess_annotations_data(annotations_data, reference_data):
    """
    Processes GTF data to extract genes, exons, introns, and intergenic regions.

    Parameters:
    - annotations_data (pd.DataFrame): GTF DataFrame, corresponding to a chromosome.
    - reference_data (dict of sequences): reference data mapping chromosome name to sequence.

    Returns:
    - dict: A dictionary containing processed data:
        - 'exons': List of exons with chromosome, strand, gene_id, gene_name, start, end, transcript_id.
        - 'introns': List of introns with chromosome, strand, gene_id, gene_name, start, end, transcript_id.
        - 'intergenic': List of intergenic regions with chromosome, start, end.
        - 'gene_annotations': List of gene annotations with chromosome, strand, gene_id, gene_name, start, end, midpoint.
    """
    # Initialize lists to store processed data
    exons_list = []
    introns_list = []
    intergenic_list = []
    gene_annotations = []

    # gtf_combined = pd.concat([annotations_data], ignore_index=True)
    gtf_combined = pd.DataFrame(annotations_data)

    # Function to extract attributes from the 'attribute' column
    def extract_attribute(attr_str, key):
        match = re.search(f'{key} "([^"]+)"', attr_str)
        return match.group(1) if match else '-' # unknown

    # Extract necessary attributes
    gtf_combined['gene_id'] = gtf_combined['attribute'].apply(lambda x: extract_attribute(x, 'gene_id'))
    gtf_combined['gene_name'] = gtf_combined['attribute'].apply(lambda x: extract_attribute(x, 'gene_name'))
    gtf_combined['transcript_id'] = gtf_combined['attribute'].apply(lambda x: extract_attribute(x, 'transcript_id'))

    # Filter to include only 'gene' features for gene annotations
    gene_df = gtf_combined[gtf_combined['feature'] == 'gene'].copy()

    # Process gene annotations
    for _, gene in gene_df.iterrows():
        gene_annotations.append({
            'chromosome': gene['seqname'],
            'strand': gene['strand'],
            'gene_id': gene['gene_id'],
            'gene_name': gene['gene_name'],
            'start': gene['start'],
            'end': gene['end'],
            'midpoint': (gene['start'] + gene['end']) / 2
        })

    # Filter to include only exons
    exon_df = gtf_combined[gtf_combined['feature'] == 'exon'].copy()

    # Sort chromosomes naturally
    chromosomes = natsorted(exon_df['seqname'].unique())

    # Process each chromosome and strand
    for chrom in chromosomes:
        print('Processing chromosome', chrom)
        chrom_exons = exon_df[exon_df['seqname'] == chrom]
        strands = ['+', '-']
        for strand in strands:
            strand_exons = chrom_exons[chrom_exons['strand'] == strand]
            if strand_exons.empty:
                continue  # Skip if no exons on this strand

            # Group exons by gene_id
            gene_groups = strand_exons.groupby('gene_id')
            gene_regions = []

            for gene_id, group in gene_groups:
                # Retrieve gene annotation for this gene_id
                gene_info = next((gene for gene in gene_annotations if gene['gene_id'] == gene_id), None)
                if gene_info is None:
                    gene_name = '-' # pseudo / unknown
                    gene_start = group['start'].min()
                    gene_end = group['end'].max()
                else:
                    gene_name = gene_info['gene_name']
                    gene_start = gene_info['start']
                    gene_end = gene_info['end']

                gene_regions.append((gene_start, gene_end, gene_id))

                # Process exons and introns within the gene
                transcripts = group['transcript_id'].unique()
                for transcript_id in transcripts:
                    transcript_exons = group[group['transcript_id'] == transcript_id].sort_values('start')
                    exons = transcript_exons[['start', 'end']].values

                    # Add exons
                    for i in range(len(exons)):
                        exon_start, exon_end = exons[i]
                        exons_list.append({
                            'chromosome': chrom,
                            'strand': strand,
                            'gene_id': gene_id,
                            'gene_name': gene_name,
                            'transcript_id': transcript_id,
                            'start': exon_start,
                            'end': exon_end
                        })

                        # Compute introns between exons
                        if i > 0:
                            prev_exon_end = exons[i - 1][1]
                            intron_start = prev_exon_end + 1
                            intron_end = exon_start - 1
                            if intron_start <= intron_end:
                                introns_list.append({
                                    'chromosome': chrom,
                                    'strand': strand,
                                    'gene_id': gene_id,
                                    'gene_name': gene_name,
                                    'transcript_id': transcript_id,
                                    'start': intron_start,
                                    'end': intron_end
                                })

            # Compute intergenic regions directly from the sorted gene_annotations
            # Get all genes on this chromosome, sorted by start
            genes_on_chrom = sorted(
                [gene for gene in gene_annotations if gene['chromosome'] == chrom and gene['strand'] == strand],
                key=lambda x: x['start']
            )

            prev_end = None
            chrom_start = 1  # Assuming chromosomes start at position 1

            for gene in genes_on_chrom:
                gene_strand = gene['strand']
                gene_start = gene['start']
                gene_end = gene['end']

                if prev_end is None:
                    # First gene on the chromosome
                    if gene_start > chrom_start:
                        intergenic_list.append({
                            'chromosome': chrom,
                            'strand': gene_strand,
                            'start': chrom_start,
                            'end': gene_start - 1
                        })
                else:
                    # Intergenic region between previous gene and current gene
                    if prev_end + 1 <= gene_start - 1:
                        intergenic_list.append({
                            'chromosome': chrom,
                            'strand': gene_strand,
                            'start': prev_end + 1,
                            'end': gene_start - 1
                        })
                prev_end = gene_end

            # handle intergenic region after the last gene (chromosome length is known)
            chrom_length = len(reference_data[chrom])
            if prev_end < chrom_length:
                intergenic_list.append({
                    'chromosome': chrom,
                    'strand': gene_strand,
                    'start': prev_end + 1,
                    'end': chrom_length
                })

    # Compile processed data into a dictionary
    processed_data = {
        'exons': exons_list,
        'introns': introns_list,
        'intergenic': intergenic_list,
        'gene_annotations': gene_annotations
    }

    return processed_data

def load_annotations(chromosome_names, reference_data):
    """
    Loads and preprocesses annotations from a list of GTF files with caching.

    For each GTF file (representing a chromosome), the function checks if a cached
    processed object exists in the 'cache' directory. If it does, it loads the cached
    data. If not, it processes the GTF file, caches the result, and then uses it.

    Parameters:
    - chromosome_names (list of str): List of chromosome names.
    - reference_data (dict of str): mapping chromosome name to sequence.

    Returns:
    - dict: Combined processed data from all chromosomes, containing:
        - 'exons'
        - 'introns'
        - 'intergenic'
        - 'gene_annotations'
    """
    cache_dir = 'cache'
    os.makedirs(cache_dir, exist_ok=True)  # Ensure cache directory exists

    combined_processed_data = {
        'exons': [],
        'introns': [],
        'intergenic': [],
        'gene_annotations': []
    }

    for chrom in chromosome_names:
        annotations_file = get_annotations_filename(chrom)
        # Derive cache file path
        base_filename = os.path.basename(annotations_file)
        cache_filename = os.path.splitext(base_filename)[0] + '.pkl'
        cache_path = os.path.join(cache_dir, cache_filename)

        if os.path.exists(cache_path):
            # Load cached processed data
            print(f"Loading cached data for {base_filename} from {cache_path}")
            with open(cache_path, 'rb') as cache_file:
                chrom_data = pickle.load(cache_file)
        else:
            # Load GTF file
            print(f"Processing annotations file: {annotations_file}")
            gtf_df = load_gtf.load_gtf(annotations_file)
            if gtf_df is None:
                print(f"Skipping {annotations_file} due to load failure.")
                continue

            # Preprocess data (assuming one chromosome per GTF file)
            chrom_data = preprocess_annotations_data(gtf_df, reference_data)

            # Save processed data to cache
            print(f"Caching processed data to {cache_path}")
            with open(cache_path, 'wb') as cache_file:
                pickle.dump(chrom_data, cache_file)

        # Combine processed data
        for key in combined_processed_data:
            combined_processed_data[key].extend(chrom_data.get(key, []))

    return combined_processed_data
