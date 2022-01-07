import pandas as pd


def load_informationcontent_from_file(path=None):
    if path is None:
        import gzip
        try:
            from importlib.resources import files
        except ImportError:
            from importlib_resources import files

        fpath = files('omamo') / "data" / "go_ic.dat"
        with gzip.open(fpath, 'rt') as fh:
            df = pd.read_csv(fh, sep='\t')
    else:
        df = pd.read_csv(path, sep='\t')
    # Create a dictionary where {GO term ID: informaiton content}
    ic = {row['t']: row['ic'] for row in df.to_records()}
    return ic
