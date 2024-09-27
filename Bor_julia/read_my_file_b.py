import pandas as pd
from itertools import combinations

def extract_name(column_path):
    """
    This will refactor the column name and only get the useful information.
    :param path:
    :return:
    """
    # Input: String of "F:\Julia\wheldone-DIFFPOP\04-19-24\OV3-whel-n3-F7.mzML"
    # Output: whel-n3-F7
    # Split the string on backslash and pick the last part
    file_name = column_path.split('\\')[-1]
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
    # Input: tsv data
    # Output: dataframe with no empty genes and nicer column names
    df = pd.read_csv(file_name, delimiter='\t')
    new_column_names = df.columns[:4].tolist() + [extract_name(col) for col in df.columns[4:]]
    df.columns = new_column_names
    df_cleaned = df.dropna(subset=['Genes'])

    select_column = ["Protein.Group"] + new_column_names[4:]

    df_cleaned = df[select_column]
    column_new_names = ["Protein"] + new_column_names[4:]

    df_cleaned.columns = column_new_names

    return df_cleaned


def dataframe_process_nsaf(file_name):
    """
    clean and rename the column names for easier downstream process
    :param file_name:
    :return:
    """
    # Read the tsv file
    # Input: tsv data
    # Output: dataframe with no empty genes and nicer column names
    combined_df = pd.read_csv(file_name, delimiter='\t')

    df_dmso = combined_df[["Protein ID"] + ["Protein Length"] + [dmso for dmso in combined_df.columns if "DMSO" in dmso]]
    df_peng = combined_df[["Protein ID"] + ["Protein Length"] + [peng for peng in combined_df.columns if "PenG" in peng]]
    df_moxi = combined_df[["Protein ID"] + ["Protein Length"] + [moxi for moxi in combined_df.columns if "Moxi" in moxi]]
    df_cipr = combined_df[["Protein ID"] + ["Protein Length"] + [cipr for cipr in combined_df.columns if "Cipro" in cipr]]


    return nsaf_method(seperate_by_result(df_dmso)), nsaf_method(seperate_by_result(df_peng)), nsaf_method(seperate_by_result(df_moxi)), nsaf_method(seperate_by_result(df_cipr))


def seperate_by_result(df):
    """
    We need spectral count for intensity, and
    :param df:
    :return:
    """
    columns = df.columns
    nsaf = ["Protein ID"] + ["Protein Length"] + [spec for spec in columns if "Total Spectral Count" in spec]
    df_nsaf = df[nsaf]

    return df_nsaf


def sort_df_columns(df, key_column="Protein"):
    """
    """
    key_col = df[key_column]

    columns_to_sort = df.columns.difference([key_column])

    sorted_columns = sorted(columns_to_sort, key=lambda x: int(x.split('-F')[1]))

    sorted_columns = [key_column] + sorted_columns

    df_sorted = df[sorted_columns]

    return df_sorted

def seperate_intensity(df):
    columns = df.columns

    df_dmso = df[["Protein"] + [dmso for dmso in columns if "DMSO" in dmso]]
    df_peng = df[["Protein"] + [peng for peng in columns if "PenG" in peng]]
    df_moxi = df[["Protein"] + [moxi for moxi in columns if "Moxi" in moxi]]
    df_cipr = df[["Protein"] + [cipr for cipr in columns if "Cipro" in cipr]]

    return sort_df_columns(df_dmso), sort_df_columns(df_peng), sort_df_columns(df_moxi), sort_df_columns(df_cipr)

def fill_na_with_half_min(df):
    """
    Since we have some of the empty cells. We fill them with half of the minimum in the row.
    :param df:
    :return:
    """
    categorical_col = df.iloc[:, 0]
    numeric_df = df.iloc[:, 1:]

    filled_numeric_df = numeric_df.apply(lambda row: row.fillna(row.min(skipna=True)/2), axis=1)
    df_filled = pd.concat([categorical_col, filled_numeric_df], axis=1)

    return df_filled


def map_column_name(column_name):
    parts = column_name.split()
    # Return as is if the column name doesn't follow the expected pattern
    if len(parts) < 2:
        return column_name

    prefix = parts[0]
    spec_type = ' '.join(parts[1:])

    if '_' not in prefix or len(prefix.split('_')) != 2:
        return column_name

    # Extract F, replicate number, and replicate set
    f, replicate_info = prefix.split('_')


    # Determine the new prefix
    if int(replicate_info) % 10 == 0:
        new_prefix = f"DMSO-n{int(replicate_info)//10}-{f}"
    else:
        new_prefix = f"whel-n{int(replicate_info)//10}-{f}"


    return f"{new_prefix}"


def nsaf_method(df):
    # Define columns to exclude from the division process
    exclude_columns = ['Protein ID', 'Protein Length']

    # Get columns that need to be divided
    columns_to_divide = [col for col in df.columns if col not in exclude_columns]

    # Strip whitespace from column names to avoid issues
    df.columns = df.columns.str.strip()

    # Make a copy of the DataFrame to avoid 'SettingWithCopyWarning'
    df = df.copy()

    # Ensure 'Protein Length' and other relevant columns are float64 to avoid division errors
    df['Protein Length'] = df['Protein Length'].astype('float64')

    # Explicitly cast all relevant columns to 'float64' in a loop
    for col in columns_to_divide:
        # Strip any potential whitespace from the column names
        col = col.strip()
        df[col] = df[col].astype('float64')

    # Perform division by 'Protein Length' for each relevant column
    for col in columns_to_divide:
        # Perform division and replace zero in 'Protein Length' with NaN to avoid division by zero
        df[col] = df[col] / df['Protein Length'].replace(0, float('nan'))

    return df




if __name__ == '__main__':
    df_intensity = fill_na_with_half_min(dataframe_process("report.pg_matrix.tsv"))
    df_dmso_nsaf, df_peng_nsaf, df_moxi_nsaf, df_cipr_nsaf = dataframe_process_nsaf("combined_protein.tsv")
    df_dmso_inten, df_peng_inten, df_moxi_inten, df_cipr_inten = seperate_intensity(df_intensity)

