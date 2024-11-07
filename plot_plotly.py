# plot_chromosomes_plotly_optimized.py

import plotly.graph_objects as go
from natsort import natsorted
import plotly.io as pio

MAX_GENE_NAME_TEXTS = 200

def plot_chromosomes_plotly(processed_data, width=1680, height=960):
    """
    Optimized Plotly function to plot chromosomes with genes, distinguishing intergenic regions and gene segments.
    Exons and introns are omitted by default to improve performance but can be toggled on interactively.

    Parameters:
    - processed_data (dict): Dictionary containing processed genomic data.
        Expected keys:
            - 'exons': List of exons with chromosome, strand, gene_name, start, end, transcript_id.
            - 'introns': List of introns with chromosome, strand, gene_name, start, end, transcript_id.
            - 'intergenic': List of intergenic regions with chromosome, strand, start, end.
            - 'gene_annotations': List of gene annotations with chromosome, strand, gene_name, start, end, midpoint.
    - width (int): Width of the figure in pixels.
    - height (int): Height of the figure in pixels.

    Output:
    - Displays an interactive Plotly plot of genes across chromosomes.
    """
    # Specify the renderer to ensure the plot opens in the browser
    pio.renderers.default = 'browser'  # e.g., 'notebook', 'iframe'

    # Extract processed data
    exons_list = processed_data.get('exons', [])
    introns_list = processed_data.get('introns', [])
    intergenic_list = processed_data.get('intergenic', [])
    gene_annotations = processed_data.get('gene_annotations', [])

    # Identify unique chromosomes and sort them naturally
    chromosomes = natsorted(list(set([
        exon['chromosome'] for exon in exons_list
    ] + [
        intron['chromosome'] for intron in introns_list
    ] + [
        intergenic['chromosome'] for intergenic in intergenic_list
    ] + [
        gene['chromosome'] for gene in gene_annotations
    ])))
    chrom_to_y = {chrom: idx for idx, chrom in enumerate(chromosomes)}

    y_tick_spacing = 3.0
    strand_offset = 0.15  # Offset for strand labels

    # Function to create segments
    def create_segments(data_list, chrom_to_y, y_tick_spacing, offset=0.3):
        segments = []
        for item in data_list:
            chrom = item['chromosome']
            strand = item['strand']
            y_pos = chrom_to_y[chrom] * y_tick_spacing
            if strand == '+':
                y_pos += offset
            else:
                y_pos -= offset
            segments.append({'x_start': item['start'], 'x_end': item['end'], 'y': y_pos})
        return segments

    # Create segments
    intergenic_segments = create_segments(intergenic_list, chrom_to_y, y_tick_spacing, offset=0.3)
    gene_segments = []
    for gene in gene_annotations:
        chrom = gene['chromosome']
        strand = gene['strand']
        y_pos = chrom_to_y[chrom] * y_tick_spacing
        if strand == '+':
            y_pos += 0.3
        else:
            y_pos -= 0.3
        gene_segments.append({
            'x_start': gene['start'],
            'x_end': gene['end'],
            'y': y_pos,
            'gene_name': gene['gene_name'],
            'midpoint': gene['midpoint'],
            'strand': strand
        })

    # Initialize figure
    fig = go.Figure()

    # Add Intergenic Regions
    if intergenic_segments:
        intergenic_x = []
        intergenic_y = []
        for seg in intergenic_segments:
            intergenic_x += [seg['x_start'], seg['x_end'], None]
            intergenic_y += [seg['y'], seg['y'], None]
        fig.add_trace(go.Scatter(
            x=intergenic_x,
            y=intergenic_y,
            mode='lines',
            line=dict(color='lightgrey', width=3),
            name='Intergenic',
            hoverinfo='none'
        ))

    # Add Genes
    if gene_segments:
        genes_x = []
        genes_y = []
        for seg in gene_segments:
            genes_x += [seg['x_start'], seg['x_end'], None]
            genes_y += [seg['y'], seg['y'], None]
        fig.add_trace(go.Scatter(
            x=genes_x,
            y=genes_y,
            mode='lines',
            line=dict(color='seagreen', width=5),
            name='Genes',
            hoverinfo='none'
        ))

    # Exons and Introns are omitted by default to improve performance
    # They can be toggled on via interactive buttons

    # Add Gene Annotations
    if gene_segments:
        annotation_x = [gene['midpoint'] for gene in gene_segments[:MAX_GENE_NAME_TEXTS]]
        annotation_y = [
            gene['y'] + strand_offset if gene['strand'] == '+' else gene['y'] - strand_offset
            for gene in gene_segments[:MAX_GENE_NAME_TEXTS]
        ]
        annotation_text = [gene['gene_name'] for gene in gene_segments[:MAX_GENE_NAME_TEXTS]]
        fig.add_trace(go.Scatter(
            x=annotation_x,
            y=annotation_y,
            mode='text',
            text=annotation_text,
            textposition='middle center',
            textfont=dict(size=6, color='darkgreen'),
            name='Gene Names',
            hoverinfo='text',
            visible=False  # Start with annotations hidden
        ))

    # Set y-axis with chromosomes and strands
    y_positions = []
    y_labels = []
    for chrom in chromosomes:
        base_y = chrom_to_y[chrom] * y_tick_spacing
        y_positions.append(base_y)
        y_labels.append(chrom)
    # Add strand labels as annotations
    strand_annotations = []
    for chrom in chromosomes:
        base_y = chrom_to_y[chrom] * y_tick_spacing
        strand_annotations.append(dict(x=0, y=base_y + strand_offset, text='+', showarrow=False, yanchor='bottom'))
        strand_annotations.append(dict(x=0, y=base_y - strand_offset, text='-', showarrow=False, yanchor='top'))
    fig.update_layout(
        yaxis=dict(
            tickmode='array',
            tickvals=y_positions,
            ticktext=y_labels,
            title="Chromosome",
            range=[min(y_positions) - y_tick_spacing, max(y_positions) + y_tick_spacing],
            zeroline=False,
            showgrid=True,
            gridwidth=0.5,
            gridcolor='lightgrey'
        ),
        annotations=strand_annotations,
        xaxis=dict(
            title="Genomic Position (bp)",
            zeroline=False,
            showgrid=True,
            gridwidth=0.5,
            gridcolor='lightgrey',
            range=[0, max(
                [item['x_end'] for item in gene_segments] +
                [item['x_end'] for item in intergenic_segments] +
                [1]  # Default if no data
            )]
        ),
        title="Genes Across Chromosomes by Strand",
        width=width,
        height=height,
        hovermode='closest'
    )

    # Add custom legend manually
    legend_items = [
        {"color": "seagreen", "label": "Genes"},
        {"color": "lightgrey", "label": "Intergenic"}
    ]
    for item in legend_items:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(size=10, color=item["color"]),
            legendgroup=item["label"],
            showlegend=True,
            name=item["label"]
        ))

    # Add Exons Trace (initially hidden)
    if exons_list:
        exon_segments = create_segments(exons_list, chrom_to_y, y_tick_spacing, offset=0.3)
        exon_x = []
        exon_y = []
        for seg in exon_segments:
            exon_x += [seg['x_start'], seg['x_end'], None]
            exon_y += [seg['y'], seg['y'], None]
        fig.add_trace(go.Scatter(
            x=exon_x,
            y=exon_y,
            mode='lines',
            line=dict(color='blue', width=10),
            name='Exons',
            hoverinfo='none',
            visible=False  # Omitted by default
        ))

    # Add Introns Trace (initially hidden)
    if introns_list:
        intron_segments = create_segments(introns_list, chrom_to_y, y_tick_spacing, offset=0.3)
        intron_x = []
        intron_y = []
        for seg in intron_segments:
            intron_x += [seg['x_start'], seg['x_end'], None]
            intron_y += [seg['y'], seg['y'], None]
        fig.add_trace(go.Scatter(
            x=intron_x,
            y=intron_y,
            mode='lines',
            line=dict(color='lightblue', width=5),
            name='Introns',
            hoverinfo='none',
            visible=False  # Omitted by default
        ))

    # Add Gene Names Trace (initially hidden)
    if gene_segments:
        fig.add_trace(go.Scatter(
            x=annotation_x,
            y=annotation_y,
            mode='text',
            text=annotation_text,
            textposition='middle center',
            textfont=dict(size=6, color='darkgreen'),
            name='Gene Names',
            hoverinfo='text',
            visible=False  # Start with annotations hidden
        ))

    # Define buttons for toggling features
    # The order of traces is important here
    # Update the visibility list accordingly
    # Trace order:
    # 0: Intergenic
    # 1: Genes
    # 2: Gene Names
    # 3: Exons (if exists)
    # 4: Introns (if exists)
    # 5+: Legend items

    # Determine the indices for exons and introns
    trace_count = len(fig.data)
    exon_trace_index = None
    intron_trace_index = None
    gene_names_trace_index = None

    for idx, trace in enumerate(fig.data):
        if trace.name == 'Exons':
            exon_trace_index = idx
        if trace.name == 'Introns':
            intron_trace_index = idx
        if trace.name == 'Gene Names':
            gene_names_trace_index = idx

    # Create buttons
    buttons = [
        dict(
            label="Show Exons",
            method="update",
            args=[
                {"visible": [True, True, False, True, False] if exon_trace_index is not None and intron_trace_index is not None else [True, True, False]},
                {"title": "Genes and Exons"}
            ]
        ),
        dict(
            label="Show Introns",
            method="update",
            args=[
                {"visible": [True, True, False, False, True] if exon_trace_index is not None and intron_trace_index is not None else [True, True, False]},
                {"title": "Genes and Introns"}
            ]
        ),
        dict(
            label="Show Gene Names",
            method="update",
            args=[
                {"visible": [True, True, True, False, False] if exon_trace_index is not None and intron_trace_index is not None else [True, True, True]},
                {"title": "Genes with Gene Names"}
            ]
        ),
        dict(
            label="Show All",
            method="update",
            args=[
                {"visible": [True, True, True, True, True] if exon_trace_index is not None and intron_trace_index is not None else [True, True, True]},
                {"title": "All Features"}
            ]
        ),
        dict(
            label="Hide All",
            method="update",
            args=[
                {"visible": [True, True, False, False, False] if exon_trace_index is not None and intron_trace_index is not None else [True, True, False]},
                {"title": "Hide All Features"}
            ]
        )
    ]

    # Add updatemenus (buttons)
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=buttons,
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.85,
                y=1.15,
                xanchor="left",
                yanchor="top"
            )
        ]
    )

    # Show the plot
    fig.show()
