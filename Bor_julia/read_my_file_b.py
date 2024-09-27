import pandas as pd
from itertools import combinations
from scipy.stats import ks_2samp
from normalization_b import *

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
    df_cleaned = df.dropna(subset=['Gene'])

    select_column = ["Gene"] + new_column_names[4:]

    df_cleaned = df[select_column]
    # column_new_names = ["Protein"] + new_column_names[4:]
    #
    # df_cleaned.columns = column_new_names

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

    df_dmso = combined_df[["Gene"] + ["Protein Length"] + [dmso for dmso in combined_df.columns if "DMSO" in dmso]]
    df_peng = combined_df[["Gene"] + ["Protein Length"] + [peng for peng in combined_df.columns if "PenG" in peng]]
    df_moxi = combined_df[["Gene"] + ["Protein Length"] + [moxi for moxi in combined_df.columns if "Moxi" in moxi]]
    df_cipr = combined_df[["Gene"] + ["Protein Length"] + [cipr for cipr in combined_df.columns if "Cipro" in cipr]]

    dmso = nsaf_method(seperate_by_result(df_dmso))


    return (new_columns(nsaf_method(seperate_by_result(df_dmso))),
            new_columns(nsaf_method(seperate_by_result(df_peng))),
            new_columns(nsaf_method(seperate_by_result(df_moxi))),
            new_columns(nsaf_method(seperate_by_result(df_cipr))))


def new_columns(df):
    columns = df.columns[1:]
    df_columns = ["-".join(name.split(" ")[0].split("_")[1:3]) for name in columns]
    df.columns = ["Gene"] + df_columns

    return df

def seperate_by_result(df):
    """
    We need spectral count for intensity, and
    :param df:
    :return:
    """
    columns = df.columns
    nsaf = ["Gene"] + ["Protein Length"] + [spec for spec in columns if "Total Spectral Count" in spec]
    df_nsaf = df[nsaf]

    return df_nsaf


def sort_df_columns(df, key_column="Gene"):
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

    df_dmso = df[["Gene"] + [dmso for dmso in columns if "DMSO" in dmso]]
    df_peng = df[["Gene"] + [peng for peng in columns if "PenG" in peng]]
    df_moxi = df[["Gene"] + [moxi for moxi in columns if "Moxi" in moxi]]
    df_cipr = df[["Gene"] + [cipr for cipr in columns if "Cipro" in cipr]]

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


def nsaf_method(df):
    exclude_columns = ['Gene', 'Protein Length']

    df.columns = df.columns.str.strip()

    columns_to_divide = [col for col in df.columns if col not in exclude_columns]

    df['Protein Length'] = df['Protein Length'].astype('float64')
    df[columns_to_divide] = df[columns_to_divide].astype('float64')

    # Step 1: Normalize by protein length
    df[columns_to_divide] = df[columns_to_divide].div(df['Protein Length'].replace(0, float('nan')), axis=0)

    # Step 2: Normalize by total abundance within each sample
    total_sums = df[columns_to_divide].sum()

    # Divide the length-normalized spectral counts by the total sum for each column (sample)
    df[columns_to_divide] = df[columns_to_divide].div(total_sums, axis=1)

    ordered_columns = ['Gene'] + columns_to_divide
    df = df[ordered_columns]

    return df


def difference_test(dmso, drug):
    print("HI")



if __name__ == '__main__':
    df_intensity = fill_na_with_half_min(dataframe_process("report.pg_matrix.tsv"))

    df_dmso_nsaf, df_peng_nsaf, df_moxi_nsaf, df_cipr_nsaf = dataframe_process_nsaf("combined_protein.tsv")

    df_dmso_inten, df_peng_inten, df_moxi_inten, df_cipr_inten = seperate_intensity(df_intensity)
    dmso_z, dmso_tic, dmso_med, dmso_quan, dmso_var = all_norm_methods(df_dmso_inten)
    peng_z, peng_tic, peng_med, peng_quan, peng_var = all_norm_methods(df_peng_inten)
    moxi_z, moxi_tic, moxi_med, moxi_quan, moxi_var = all_norm_methods(df_moxi_inten)
    cipr_z, cipr_tic, cipr_med, cipr_quan, cipr_var = all_norm_methods(df_cipr_inten)

    plot_gene_intensity_DP(df_dmso_nsaf, df_peng_nsaf,  "BB_0732", "NSAF")
    plot_gene_intensity_CM(df_dmso_nsaf, df_cipr_nsaf,  "gyrB", "NSAF")
    plot_gene_intensity_DM(df_dmso_nsaf, df_moxi_nsaf,  "gyrB", "NSAF")
    plot_gene_intensity_CM(df_dmso_nsaf, df_cipr_nsaf,  "gyrA", "NSAF")
    plot_gene_intensity_DM(df_dmso_nsaf, df_moxi_nsaf,  "gyrA", "NSAF")

    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplA", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplB", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplC", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplD", "NSAF")
    # plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsA", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsB", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsC", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsD", "NSAF")

