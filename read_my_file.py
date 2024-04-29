import pandas as pd
from itertools import combinations

def extract_name(path):
    """
    This will refactor the column name and only get the useful information.
    :param path:
    :return:
    """
    # Split the string on backslash and pick the last part
    file_name = path.split('\\')[-1]
    # Split the file name on '-' and extract the necessary parts
    parts = file_name.split('-')
    with_mzml = '-'.join(parts[1:4])
    return with_mzml.split('.')[0]

def dataframe_process(file_name):
    """
    clean and rename the column names for easier downstream process
    :param file_name:
    :return:
    """
    # Read the tsv file
    df = pd.read_csv(file_name, delimiter='\t')
    new_column_names = df.columns[:5].tolist() + [extract_name(col) for col in df.columns[5:]]
    df.columns = new_column_names
    df_cleaned = df.dropna(subset=['Genes'])

    return df_cleaned


def parse_col_name(col):
    """
    This will parse the column name that can tell us whether it is whel or DSMO, F1 to F9 and n1 to n3
    :param col:
    :return:
    """
    prefix, run, f = col.split('-')
    return prefix, run, f


def column_group(df):
    """
    get the dataframe of the group columns for pearson correlation
    :param df:
    :return:
    """
    groups = {}
    for col in df.columns[5:]:
        prefix, run, f_number = parse_col_name(col)
        key = (prefix, f_number)
        if key not in groups:
            groups[key] = []
        groups[key].append(col)

    return groups


def pearson_correlations(df, group):
    """
    Get the subset correlation matrix
    :param df:
    :param group:
    :return:
    """
    # Extract only the relevant columns for the group
    subset = df[group]
    # Compute the correlation matrix
    return subset.corr()


def dataframe_work(df, column_dict):
    """
    get the dictionary of all the pearson correlations.
    :param df:
    :param column_dict:
    :return:
    """

    results = {}

    for key, group_columns in column_dict.items():
            results[key] = pearson_correlations(df, group_columns)
            print(f"Correlations for {key}:")
            print(results[key])

    return results

if __name__ == '__main__':
    dd = dataframe_process("report.pg_matrix 3.tsv")
    gg = column_group(dd)
    aa = dataframe_work(dd, gg)
