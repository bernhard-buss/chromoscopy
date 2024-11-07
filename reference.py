from Bio import SeqIO
import os
import pickle

reference_base_name = 'reference/Homo_sapiens.GRCh38.dna.chromosome.'
cache_dir = 'cache'

def get_reference_chromosome_filename(chrom):
    return f"{reference_base_name}{chrom}.fa"

def load_chromosome_sequence(chrom):
    filename = get_reference_chromosome_filename(chrom)
    base_filename = os.path.basename(filename)
    cache_filename = os.path.splitext(base_filename)[0] + '.pkl'
    cache_path = os.path.join(cache_dir, cache_filename)

    if os.path.exists(cache_path):
        # Load cached reference data
        print(f"Loading cached data for {base_filename} from {cache_path}")
        with open(cache_path, 'rb') as cache_file:
            seq = pickle.load(cache_file)
    else:
        print('Reading reference sequence for chrom', chrom)
        seq = SeqIO.read(filename, "fasta").seq  # Load only the sequence

        # Save processed data to cache
        print(f"Caching processed data to {cache_path}")
        with open(cache_path, 'wb') as cache_file:
            pickle.dump(seq, cache_file)

    return seq

def load_reference_data(chromosome_names):
    reference_sequences = [load_chromosome_sequence(chrom) for chrom in chromosome_names]
    reference_data = dict()
    for i in range(len(reference_sequences)):
        reference_data[chromosome_names[i]] = reference_sequences[i]

    return reference_data
