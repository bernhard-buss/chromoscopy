# plot_chromosomes.py

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter, AutoLocator
from matplotlib.widgets import CheckButtons, TextBox
import matplotlib.patches as mpatches
from natsort import natsorted
from collections import defaultdict

MAX_GENE_NAME_TEXTS = 200

def plot_chromosomes(reference_data, annotations_data, dpi=100, width_pixels=1680, height_pixels=960, y_tick_spacing=3.0):
    """
    Plots chromosomes with genes, distinguishing exons, introns, intergenic regions, and gene segments.
    Includes gene name annotations that dynamically show/hide based on zoom level and visible strands.
    Adds '+' and '-' annotations on the y-axis, with chromosomes as main ticks and strands as sub-ticks.
    Provides interactive widgets to toggle visibility of plot features.
    Ensures the plot cannot be panned below x = 0.

    Parameters:
    - reference_data (dict): maps chromosome names to sequence
    - annotations_data (dict): Dictionary containing processed genomic data from preprocess_gtf_data.
        Expected keys:
            - 'exons': List of exons with chromosome, strand, gene_id, gene_name, start, end, transcript_id.
            - 'introns': List of introns with chromosome, strand, gene_id, gene_name, start, end, transcript_id.
            - 'intergenic': List of intergenic regions with chromosome, strand, start, end.
            - 'gene_annotations': List of gene annotations with chromosome, strand, gene_id, gene_name, start, end, midpoint.
    - dpi (int): Dots per inch for the figure resolution.
    - width_pixels (int): Desired width of the figure in pixels.
    - height_pixels (int): Desired height of the figure in pixels.
    - y_tick_spacing (float): Spacing between chromosomes on the y-axis.

    Output:
    - Displays an interactive matplotlib plot of genes across chromosomes.
    """
    # Extract processed data
    exons_list = annotations_data.get('exons', [])
    introns_list = annotations_data.get('introns', [])
    intergenic_list = annotations_data.get('intergenic', [])
    gene_annotations = annotations_data.get('gene_annotations', [])

    # Preprocess exons and introns to get transcripts
    gene_to_transcripts = defaultdict(list)
    transcript_dict = defaultdict(lambda: {'start': float('inf'), 'end': float('-inf'), 'exons': [], 'introns': []})

    # Group exons by gene_id and transcript_id
    for exon in exons_list:
        gene_id = exon['gene_id']
        transcript_id = exon.get('transcript_id', 'unknown')
        key = (gene_id, transcript_id)
        transcript = transcript_dict[key]
        transcript['start'] = min(transcript['start'], exon['start'])
        transcript['end'] = max(transcript['end'], exon['end'])
        transcript['exons'].append((exon['start'], exon['end']))

    # Group introns by gene_id and transcript_id
    for intron in introns_list:
        gene_id = intron['gene_id']
        transcript_id = intron.get('transcript_id', 'unknown')
        key = (gene_id, transcript_id)
        transcript = transcript_dict[key]
        transcript['introns'].append((intron['start'], intron['end']))

    # Populate gene_to_transcripts
    for (gene_id, transcript_id), transcript in transcript_dict.items():
        gene_to_transcripts[gene_id].append({
            'transcript_id': transcript_id,
            'start': transcript['start'],
            'end': transcript['end'],
            'exons': transcript['exons'],
            'introns': transcript['introns']
        })

    # Identify unique chromosomes and sort them naturally
    chromosomes = natsorted(list(set([exon['chromosome'] for exon in exons_list])))
    chrom_to_y = {chrom: idx for idx, chrom in enumerate(chromosomes)}

    def strand_y_pos(chrom, strand):
        y_pos = chrom_to_y[chrom] * y_tick_spacing
        if strand == '+':
            y_pos += y_tick_spacing * 0.1  # Adjusting position for '+' strand
        else:
            y_pos -= y_tick_spacing * 0.1  # Adjusting position for '-' strand
        return y_pos

    # Initialize lists to store plotting segments
    exon_segments = []
    intron_segments = []
    intergenic_segments = []
    gene_segments = []

    # Organize exons into segments for LineCollection
    def compute_exon_segments():
        if exon_segments:
            return
        for exon in exons_list:
            y_pos = strand_y_pos(exon['chromosome'], exon['strand'])
            exon_segments.append([(exon['start'], y_pos), (exon['end'], y_pos)])

    # Organize introns into segments for LineCollection
    def compute_intron_segments():
        if intron_segments:
            return
        for intron in introns_list:
            y_pos = strand_y_pos(intron['chromosome'], intron['strand'])
            intron_segments.append([(intron['start'], y_pos), (intron['end'], y_pos)])

    # Compute exon and intron segments
    compute_exon_segments()
    compute_intron_segments()

    # Organize intergenic regions into segments
    for intergenic in intergenic_list:
        y_pos = strand_y_pos(intergenic['chromosome'], intergenic['strand'])
        intergenic_segments.append([(intergenic['start'], y_pos), (intergenic['end'], y_pos)])

    # Organize gene_annotations into gene_segments
    for gene in gene_annotations:
        y_pos = strand_y_pos(gene['chromosome'], gene['strand'])
        gene_segments.append([(gene['start'], y_pos), (gene['end'], y_pos)])

    # Set up the figure
    width_in_inches = width_pixels / dpi
    height_in_inches = height_pixels / dpi
    fig_height = max(2, len(chromosomes) * y_tick_spacing)

    fig, ax = plt.subplots(figsize=(width_in_inches, min(height_in_inches, fig_height)), dpi=dpi)

    # Create LineCollections for intergenic regions
    if intergenic_segments:
        intergenic_lc = LineCollection(intergenic_segments, colors='lightgrey', linewidths=3, label='Intergenic', visible=True)
        ax.add_collection(intergenic_lc)
    else:
        intergenic_lc = None

    # Create LineCollections for genes, enabling picking
    if gene_segments:
        gene_lc = LineCollection(gene_segments, colors='seagreen', linewidths=5, label='Genes', visible=True, picker=True)
        ax.add_collection(gene_lc)
        # Attach gene information to the LineCollection for easy access in the event handler
        gene_lc.genes = gene_annotations
    else:
        gene_lc = None

    # Create LineCollections for introns and exons, initially hidden
    intron_lc = None
    def add_intron_lines():
        nonlocal intron_lc
        if intron_lc:
            return
        if intron_segments:
            intron_lc = LineCollection(intron_segments, colors='lightblue', linewidths=5, label='Introns', visible=False)
            ax.add_collection(intron_lc)
        else:
            intron_lc = None

    exon_lc = None
    def add_exon_lines():
        nonlocal exon_lc
        if exon_lc:
            return
        if exon_segments:
            exon_lc = LineCollection(exon_segments, colors='blue', linewidths=10, label='Exons', visible=False)
            ax.add_collection(exon_lc)
        else:
            exon_lc = None

    # Assign y positions to chromosomes and strands
    y_positions = []
    y_labels = []
    minor_y_positions = []
    minor_y_labels = []
    for chrom in chromosomes:
        base_y = chrom_to_y[chrom] * y_tick_spacing
        y_positions.append(base_y)
        y_labels.append(chrom)
        # '+' and '-' strands
        minor_y_positions.extend([strand_y_pos(chrom, '+'), strand_y_pos(chrom, '-')])
        minor_y_labels.extend(['+', '-'])

    # Determine x-axis limits ensuring x_min >= 0
    all_x_positions = []
    for segment in gene_segments + intergenic_segments:
        all_x_positions.extend([pos[0] for pos in segment] + [pos[1] for pos in segment])
    if all_x_positions:
        start_min = max(min(all_x_positions), 0)
        end_max = max(all_x_positions)
        ax.set_xlim(start_min, end_max)
    else:
        ax.set_xlim(0, 1)  # Default limits if no data

    # ============================
    # Customize Plot Aesthetics
    # ============================
    ax.set_title('Genes Across Chromosomes by Strand', fontsize=14)
    ax.set_xlabel('Genomic Position (bp)', fontsize=12)
    ax.set_ylabel('Chromosome and Strand', fontsize=12)

    def abbreviate_number(x, pos):
        if x >= 1e9:
            return f'{x*1e-9:.1f}G'
        elif x >= 1e6:
            return f'{x*1e-6:.1f}M'
        #elif x >= 1e3:
        #    return f'{x*1e-3:.3f}k'
        else:
            return f'{x:.0f}'
    # Set x-axis formatter and locator
    ax.xaxis.set_major_formatter(FuncFormatter(abbreviate_number))
    ax.xaxis.set_major_locator(AutoLocator())

    # Set y-ticks
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=10)
    # Adjust the padding between y-axis tick labels and the axis
    ax.tick_params(axis='y', which='major', pad=15)

    # Set minor y-ticks for strands
    ax.set_yticks(minor_y_positions, minor=True)
    ax.set_yticklabels(minor_y_labels, minor=True, fontsize=8)

    # Set y-axis limits
    ax.set_ylim(min(y_positions) - y_tick_spacing * 0.5, max(y_positions) + y_tick_spacing * 0.5)

    # Add grid lines for x-axis
    ax.grid(axis='x', linestyle='--', alpha=0.7)

    # Create CheckButtons for toggling features
    rax = plt.axes([0.85, 0.4, 0.1, 0.2])  # [left, bottom, width, height]
    labels = ['Genes', 'Exons', 'Introns', 'Intergenic', 'Gene Names']
    visibility = [True, False, False, True, False]  # 'Genes' checked by default; 'Exons' and 'Introns' not checked
    check = CheckButtons(rax, labels, visibility)

    # Define the position and size of the search textbox (left, bottom, width, height)
    search_text_box_ax = plt.axes([0.85, 0.68, 0.1, 0.05])  # aligned with checkboxes on the right
    search_text_box_ax.text(0, 1.4, 'Search Gene:', transform=search_text_box_ax.transAxes,
                 fontsize=10, verticalalignment='center', horizontalalignment='left')

    # Create the TextBox widget
    search_gene_text_box = TextBox(search_text_box_ax, '', initial='', color='lightyellow', hovercolor='white')

    # Initialize a list to keep track of currently visible Text objects (gene names)
    current_visible_texts = []
    # Initialize a dictionary to keep track of plotted transcripts for each gene_id
    plotted_transcripts = {}

    # Function to update gene names based on current view
    def update_gene_names(x_min, x_max, y_min, y_max):
        nonlocal current_visible_texts
        # Remove existing annotations
        for text in current_visible_texts:
            text.remove()
        current_visible_texts = []

        if not check.get_status()[4]:  # 'Gene Names' checkbox is unchecked
            return

        # Iterate over gene_annotations and add annotations within current view
        current_visible_genes = []
        for gene in gene_annotations:
            if (gene['start'] <= x_max and gene['end'] >= x_min and
                y_min <= strand_y_pos(gene['chromosome'], gene['strand']) <= y_max):
                current_visible_genes.append(gene)
        if len(current_visible_genes) < MAX_GENE_NAME_TEXTS:
            for gene in current_visible_genes:
                # Assign y position based on strand
                y_pos = strand_y_pos(gene['chromosome'], gene['strand'])
                if gene['strand'] == '+':
                    y_text = y_pos + y_tick_spacing * 0.02
                    va = 'bottom'
                else:
                    y_text = y_pos - y_tick_spacing * 0.02
                    va = 'top'

                # Add text annotation using gene_name
                text = ax.text(
                    gene['midpoint'], y_text, gene['gene_name'],
                    ha='center', va=va, fontsize=8, rotation=45, color='darkgreen'
                )
                current_visible_texts.append(text)

    # Initialize a list to keep track of currently visible sequence Text objects
    current_visible_sequences = []

    def update_sequence_text(x_min, x_max, y_min, y_max):
        nonlocal current_visible_sequences

        # Remove existing sequence texts
        for text in current_visible_sequences:
            text.remove()
        current_visible_sequences.clear()

        # Only render sequence if x-range is less than 200 bp
        if (x_max - x_min) > 200:
            return

        # Adjust font size to match the line height of the intergenic region
        # Assuming that the intergenic region is represented by 'intergenic_lc'
        # We'll estimate the line height based on y_tick_spacing
        line_height = y_tick_spacing * 0.05  # Adjust as needed
        font_size = line_height * 72  # Convert line height to points (1 point = 1/72 inches)

        # Iterate over chromosomes visible in the current y-range
        for chrom_index, _ in enumerate(y_positions):
            chromosome = y_labels[chrom_index]  # Get chromosome name
            strand = '+'  # we consider only the '+' strand
            y_pos_strand = strand_y_pos(chromosome, strand)
            if y_min <= y_pos_strand <= y_max:
                # Get the sequence for the current x-range
                sequence = reference_data.get(chromosome, '')  # Get sequence for chromosome
                if not sequence:
                    continue  # No sequence data for this chromosome

                # Extract the sequence in the current x-range
                seq_start = int(max(0, x_min))
                seq_end = int(min(len(sequence), x_max))
                seq_fragment = sequence[seq_start:seq_end]

                if not seq_fragment:
                    continue  # No sequence data in this range

                # Render the sequence text
                # Position the text at the y_pos_strand corresponding to the chromosome and strand
                for i in range(len(seq_fragment)):
                    text = ax.text(
                        seq_start+i+1, y_pos_strand, seq_fragment[i],
                        fontsize=font_size, ha='center', va='center',
                        fontfamily='monospace', color='black'
                    )
                    current_visible_sequences.append(text)


    # Callback function for CheckButtons
    def toggle_features(label):
        if label == 'Genes' and gene_lc:
            gene_lc.set_visible(not gene_lc.get_visible())
        elif label == 'Exons':
            compute_exon_segments()
            add_exon_lines()
            if exon_lc:
                exon_lc.set_visible(not exon_lc.get_visible())
        elif label == 'Introns':
            compute_intron_segments()
            add_intron_lines()
            if intron_lc:
                intron_lc.set_visible(not intron_lc.get_visible())
        elif label == 'Intergenic' and intergenic_lc:
            intergenic_lc.set_visible(not intergenic_lc.get_visible())
        elif label == 'Gene Names':
            x_min, x_max = ax.get_xlim()
            y_min, y_max = ax.get_ylim()
            if check.get_status()[4]:
                update_gene_names(x_min, x_max, y_min, y_max)
            else:
                # Remove all annotations
                for text in current_visible_texts:
                    text.remove()
                current_visible_texts.clear()
        plt.draw()

    # Connect the callback to CheckButtons
    check.on_clicked(toggle_features)

    def on_text_submit(text):
        gene_name = text.strip()
        # Search for the gene in your annotations
        matching_genes = [gene for gene in gene_annotations if gene['gene_name'] == gene_name]
        if not matching_genes:
            print(f"Gene '{gene_name}' not found.")
            return
        # Assuming gene names are unique, take the first match
        gene = matching_genes[0]
        # Get gene position and chromosome
        gene_start = gene['start']
        gene_end = gene['end']
        chromosome = gene['chromosome']
        strand = gene['strand']
        # Calculate new x-limits with some padding
        padding = (gene_end - gene_start) * 0.1  # 10% padding
        x_min_new = max(0, gene_start - padding)
        x_max_new = gene_end + padding
        # Calculate new y-limits to center on the chromosome and strand
        y_pos = strand_y_pos(chromosome, strand)
        y_tick_spacing = y_positions[1] - y_positions[0]  # Assuming uniform spacing
        y_min_new = y_pos - y_tick_spacing * 1.5
        y_max_new = y_pos + y_tick_spacing * 1.5
        # Set the new limits
        ax.set_xlim(x_min_new, x_max_new)
        ax.set_ylim(y_min_new, y_max_new)
        # Redraw the plot
        plt.draw()

    search_gene_text_box.on_submit(on_text_submit)

    # Function to handle axis limit changes (zoom/pan)
    def on_limits_changed(event_ax):
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        # Prevent x_min from going below 0
        if x_min < 0:
            ax.set_xlim(0, x_max - x_min)
            return  # Exit to prevent further processing

        # Enforce minimum y-range (e.g., 2 chromosomes)
        min_chromosomes = 2
        min_y_range = y_tick_spacing * (min_chromosomes - 1)
        current_y_range = y_max - y_min

        if current_y_range < min_y_range:
            # Calculate new y-limits centered around the current center
            y_center = (y_max + y_min) / 2
            y_min_new = y_center - min_y_range / 2
            y_max_new = y_center + min_y_range / 2

            # Ensure y-limits are within the overall data range
            y_min_total = min(y_positions) - y_tick_spacing * 0.5
            y_max_total = max(y_positions) + y_tick_spacing * 0.5

            if y_min_new < y_min_total:
                y_min_new = y_min_total
                y_max_new = y_min_new + min_y_range
            if y_max_new > y_max_total:
                y_max_new = y_max_total
                y_min_new = y_max_new - min_y_range

            ax.set_ylim(y_min_new, y_max_new)
            return  # Exit to prevent further processing

        # Enforce minimum x-range (e.g., 100 bp)
        min_x_range = 100  # Minimum x-range in bp
        current_x_range = x_max - x_min

        if current_x_range < min_x_range:
            # Calculate new x-limits centered around the current center
            x_center = (x_max + x_min) / 2
            x_min_new = x_center - min_x_range / 2
            x_max_new = x_center + min_x_range / 2

            # Ensure x-limits are within the overall data range
            x_min_total = 0  # Assuming genomic positions start at 0
            if x_min_new < x_min_total:
                x_min_new = x_min_total
                x_max_new = x_min_new + min_x_range

            ax.set_xlim(x_min_new, x_max_new)
            return  # Exit to prevent further processing

        # Update gene names based on new limits
        update_gene_names(x_min, x_max, y_min, y_max)

        # Update sequence text based on new limits
        update_sequence_text(x_min, x_max, y_min, y_max)

    # Connect the event handler to both xlim_changed and ylim_changed events
    ax.callbacks.connect('xlim_changed', on_limits_changed)
    ax.callbacks.connect('ylim_changed', on_limits_changed)

    def on_pick(event):
        """
        Event handler for pick events. Detects if a gene line was clicked and toggles its transcripts.
        """
        # Check if the artist picked is our gene LineCollection
        if event.artist == gene_lc:
            # event.ind contains the indices of the picked lines
            for ind in event.ind:
                gene = gene_lc.genes[ind]
                gene_id = gene['gene_id']
                gene_name = gene['gene_name']
                strand = gene['strand']
                chromosome = gene['chromosome']

                if gene_id in plotted_transcripts:
                    # Transcripts are already plotted; remove them
                    for lc in plotted_transcripts[gene_id]:
                        lc.remove()
                    del plotted_transcripts[gene_id]
                else:
                    # Plot transcripts
                    transcripts = gene_to_transcripts.get(gene_id, [])
                    if not transcripts:
                        print(f"No transcripts found for {gene_name} (ID: {gene_id})")
                        continue

                    # Determine base y position
                    y_pos = strand_y_pos(chromosome, strand)
                    offset_step = y_tick_spacing * 0.01  # Minimal offset for closeness

                    # List to store LineCollections for this gene's transcripts
                    transcript_lcs = []
                    transcripts_exons_segments = []
                    transcripts_introns_segments = []

                    for i, transcript in enumerate(transcripts):
                        # Calculate y offset
                        if strand == '+':
                            y_offset = y_pos + (i + 1) * offset_step
                        else:
                            y_offset = y_pos - (i + 1) * offset_step
                        # Plot exons
                        for exon_start, exon_end in sorted(transcript['exons']):
                            transcripts_exons_segments.append([(exon_start, y_offset), (exon_end, y_offset)])
                        # Plot introns
                        for intron_start, intron_end in sorted(transcript['introns']):
                            transcripts_introns_segments.append([(intron_start, y_offset), (intron_end, y_offset)])

                    transcripts_exons_lc = LineCollection(transcripts_exons_segments, colors='blue', linewidths=5, label=f"{gene_name} Transcripts Exons", visible=True)
                    ax.add_collection(transcripts_exons_lc)
                    transcripts_introns_lc = LineCollection(transcripts_introns_segments, colors='lightblue', linewidths=3, label=f"{gene_name} Transcripts Introns", visible=True)
                    ax.add_collection(transcripts_introns_lc)
                    # Store the transcript LineCollection to allow removal later
                    transcript_lcs.append(transcripts_exons_lc)
                    transcript_lcs.append(transcripts_introns_lc)
                    # Store the transcript LineCollections
                    plotted_transcripts[gene_id] = transcript_lcs
                    print(f"Plotted {len(transcripts)} transcripts for {gene_name} (ID: {gene_id})")

                plt.draw()

    # Connect the pick event to the handler
    fig.canvas.mpl_connect('pick_event', on_pick)

    # Initially update gene names based on current view and checkbox state
    x_min_init, x_max_init = ax.get_xlim()
    y_min_init, y_max_init = ax.get_ylim()
    update_gene_names(x_min_init, x_max_init, y_min_init, y_max_init)

    # Add custom legend
    gene_patch = mpatches.Patch(color='seagreen', label='Genes')
    exon_patch = mpatches.Patch(color='blue', label='Exons')
    intron_patch = mpatches.Patch(color='lightblue', label='Introns')
    intergenic_patch = mpatches.Patch(color='lightgrey', label='Intergenic')
    ax.legend(handles=[gene_patch, exon_patch, intron_patch, intergenic_patch], loc='upper right')

    # Adjust layout to accommodate CheckButtons
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space on the right for widgets

    plt.show()
