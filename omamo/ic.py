import pandas as pd


def load_informationcontent_from_file(path):
    df = pd.read_csv(path, sep='\t')
    # Create a dictionary where {GO term ID: informaiton content}
    ic = {row['t']: row['ic'] for row in df.to_records()}
    return ic



