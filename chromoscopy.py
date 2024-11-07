import os
import pandas as pd
import annotations
import reference
import plot_interactive
import plot_genome_tracks
import plot_plotly

# ============================
# Configure Matplotlib Backend
# ============================
import matplotlib
# matplotlib.use('TkAgg')  # Use 'TkAgg' for interactive features (default: 'MacOSX')

# ============================
# Setup
# ============================
# os.chdir('chromoscopy') # if using a different working directory

# List of chromosome names
chromosome_names = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
#chromosome_names = ['Y']

# ============================
# Preprocess FASTA reference
# ============================
print('Processing reference data...')
reference_data = reference.load_reference_data(chromosome_names)

# ============================
# Preprocess GTF annotations
# ============================
print('Processing annotations data...')
annotations_data = annotations.load_annotations(chromosome_names, reference_data)

# ============================
# Generate Chromosome Plot
# ============================
print('Generating chromosome plot...')
#plot_genome_tracks.plot_genome_tracks(annotations_data)
plot_interactive.plot_chromosomes(reference_data, annotations_data)
#plot_plotly.plot_chromosomes_plotly(annotations_data)
