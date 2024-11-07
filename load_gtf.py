import pandas as pd

def load_gtf(gtf_filename):
    """
    Loads a GTF file into a Pandas DataFrame.

    Parameters:
    - gtf_filename (str): Path to the GTF file.

    Returns:
    - pd.DataFrame: DataFrame containing the GTF data.
    """

    gtf_columns = [
        'seqname', 'source', 'feature', 'start', 'end',
        'score', 'strand', 'frame', 'attribute'
    ]
    try:
        gtf_df = pd.read_csv(
            gtf_filename,
            sep='\t',
            comment='#',
            header=None,
            names=gtf_columns,
            dtype={'seqname': str}
        )
        return gtf_df
    except FileNotFoundError:
        print(f"File {gtf_filename} not found.")
        return None
