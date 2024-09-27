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


    return seperate_by_result(df_dmso), seperate_by_result(df_peng), seperate_by_result(df_moxi), seperate_by_result(df_cipr)


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


def spearman_correlations(df, group):
    """
    Get the subset correlation matrix
    :param df:
    :param group:
    :return:
    """
    # Extract only the relevant columns for the group
    subset = df[group]
    # Compute the correlation matrix
    return subset.corr(method='spearman')


def dataframe_work(df, column_dict):
    """
    get the dictionary of all the pearson correlations.
    :param df:
    :param column_dict:
    :return:
    """

    results = {}
    # Get the correlation within each fraction
    for key, group_columns in column_dict.items():
            results[key] = pearson_correlations(df, group_columns)
            print(f"Correlations for {key}:")
            print(results[key])

    return results


def dataframe_work_spearman(df, column_dict):
    """
    get the dictionary of all the spearman correlations.
    :param df:
    :param column_dict:
    :return:
    """

    results = {}
    # Get the correlation within each fraction
    for key, group_columns in column_dict.items():
            results[key] = spearman_correlations(df, group_columns)
            print(f"Correlations for {key}:")
            print(results[key])

    return results


def merge_new_run(df, df_new):
    """
    Merge the two dataframes and replacing the new runs to the old runs on the same gene.
    :param df:
    :param df_new:
    :return:
    """
    # Input: 2 dataframes old and new
    # Output: a merged dataframe with updated runs
    col_n = list(df_new.columns)[5:-1]
    merged_df = df.merge(df_new, on='Genes', how='left', suffixes=('', '_new'))
    merged_df.drop([merged_df.columns[i] + "_new" for i in [0, 1, 2, 4]], axis=1, inplace=True)
    merged_df.drop(merged_df.columns[-1], axis=1, inplace = True)
    merged_df.columns = [rename_column(col) for col in merged_df.columns]

    # for col in col_n:
    #     merged_df[col] = merged_df[col + "_new"]
    # merged_df.drop([col + '_new' for col in col_n], axis=1, inplace=True)

    return merged_df


def rename_column(col_name):
    """
    Change the column names from n(x)_new to r(x)
    :param col_name:
    :return:
    """
    if col_name.endswith('_new'):
        new_name = col_name.replace('_new', '').replace('n', 'r')
        return new_name
    return col_name


def transform_correlation_dict(corr_dict):
    """
    Get the correlation dictionary into a dataframe for output as excel or tsv.
    :param corr_dict:
    :return:
    """
    results = []

    # Iterate through the dictionary and get all the combinations in the set.
    for (condition, fraction), df in corr_dict.items():
        correlations = {'Fraction': f'{condition}-{fraction}'}
        processed_pairs = set()
        pairs = list(combinations(df.columns, 2))
        for (col1, col2) in pairs:
            pair_name = f"{col1.split('-')[1]}+{col2.split('-')[1]}"
            reverse_pair_name = f"{col2.split('-')[1]}+{col1.split('-')[1]}"
            if pair_name not in processed_pairs and reverse_pair_name not in processed_pairs:
                correlations[pair_name] = df.loc[col1, col2]
                processed_pairs.add(pair_name)
        results.append(correlations)
    final_df = pd.DataFrame(results)

    #  Still need manual adjustments.
    return final_df


def separate_dataframe(df):
    """
    THis will seperate the dataframe into two dataframes.
    :param df:
    :return:
    """
    genes_col = df['Genes']

    dmso_cols = [col for col in df.columns if col.startswith('DMSO-')]
    whel_cols = [col for col in df.columns if col.startswith('whel-')]

    df_dmso = pd.DataFrame({'Genes': genes_col})
    df_dmso = pd.concat([df_dmso, df[dmso_cols]], axis=1)

    df_whel = pd.DataFrame({'Genes': genes_col})
    df_whel = pd.concat([df_whel, df[whel_cols]], axis=1)

    return df_dmso, df_whel


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


def fill_na_with_zero(df):
    """
    Since we have some of the empty cells. We fill them with half of the minimum in the row.
    :param df:
    :return:
    """
    categorical_col = df.iloc[:, 0]
    numeric_df = df.iloc[:, 1:]

    filled_numeric_df = numeric_df.apply(lambda row: row.fillna(0), axis=1)
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


def rename_dataframe_columns(df):
    new_column_names = {col: map_column_name(col) for col in df.columns}
    # Create a copy of the DataFrame and rename the columns
    df_copy = df.copy()
    df_copy.rename(columns=new_column_names, inplace=True)

    return df_copy


# # Read all the files needed
# ccp = pd.read_excel("ccp.xlsx", "Atlas")
# ccp_genes = list(set(ccp["Genes"]))
# df_new = pd.read_csv("1-6_MBR.tsv", delimiter= "\t")
# df_new_filtered = df_new[~df_new['Protein'].str.startswith('tr|') & ~df_new['Protein'].str.contains('contam')]
# intensity_pattern = r'F\d+_\d+ Intensity'
# columns_to_subset_intensity = df_new_filtered.filter(regex=intensity_pattern).columns
# columns_to_subset_intensity = ["Genes"] + list(columns_to_subset_intensity)
# df_intensity_intermediate = df_new_filtered[columns_to_subset_intensity]
# intensity_df = rename_dataframe_columns(df_intensity_intermediate)
#
# df_dmso, df_whel = separate_dataframe(intensity_df)
# filled_dmso = fill_na_with_half_min(df_dmso).dropna()
# filled_whel = fill_na_with_half_min(df_whel).dropna()
#
# spec_count_df = pd.read_csv("1-6_MBR.tsv", delimiter= "\t")
# spec_count_df_filtered = spec_count_df[~spec_count_df['Protein'].str.startswith('tr|') & ~spec_count_df['Protein'].str.contains('contam')]
# spec_pattern = r'F\d+_\d+ Total Spectral Count'
# columns_to_subset = spec_count_df_filtered.filter(regex=spec_pattern).columns
# columns_to_subset = ["Genes"] + ["Protein Length"]+ list(columns_to_subset)
#
# # Rename the columns
# sc_df_intermediate = spec_count_df_filtered[columns_to_subset]
# spec_df = rename_dataframe_columns(sc_df_intermediate)
#
# # Get all the genes for each run and Union them
# dmso_gene = set(filled_dmso["Genes"])
# whel_gene = set(filled_whel["Genes"])
# inter_genes = list(dmso_gene & whel_gene)
#
# # Subset shared genes
# dmso_shared = filled_dmso[filled_dmso["Genes"].isin(inter_genes)]
# whel_shared = filled_whel[filled_whel["Genes"].isin(inter_genes)]

if __name__ == '__main__':
    df_intensity = dataframe_process("report.pg_matrix.tsv")
    df_dmso_nsaf, df_peng_nsaf, df_moxi_nsaf, df_cipr_nsaf = dataframe_process_nsaf("combined_protein.tsv")


