import pandas as pd

from read_my_file import *
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import numpy as np
from sklearn.preprocessing import QuantileTransformer
from pathway_construction import *
import seaborn as sns

# def whel_dmso_subset(df):
#     """
#     Subset the dataframe into gene + whel-xxx and gene + dmso-xxx
#     :param df:
#     :return:
#     """
#     gene_column = [col for col in df.columns if col == 'Genes']
#     whel_columns = [col for col in df.columns if col.startswith('whel-')]
#     dmso_columns = [col for col in df.columns if col.startswith('DMSO-')]
#
#     gene_whel_data = df[gene_column + whel_columns]
#     gene_dmso_data = df[gene_column + dmso_columns]
#
#     dmso_runs = subset_by_run(gene_dmso_data)
#     whel_runs = subset_by_run(gene_whel_data)
#
#     return whel_runs, dmso_runs


def subset_by_run(df):
    """
    Subset my dataframe by run. This looks at the run number and subset each of them into runs.
    :param df:
    :return:
    """
    genes_column = df[['Genes']]

    x_values = set(col.split('-')[1][1:] for col in df.columns if '-' in col and col != 'Genes')

    subsets = {}

    for x in x_values:
        # Since there are two types of runs. We need r and n individually.
        columns_with_x = [col for col in df.columns if f'-r{x}-' in col or f'-n{x}-' in col]
        columns_with_x_sorted = sorted(columns_with_x, key=lambda col: int(col.split('-F')[-1]))
        subset = pd.concat([genes_column, df[columns_with_x_sorted]], axis=1)
        subsets[x] = subset
    return subsets

# def graphing_methods(dataframe, gene_name,treatment_type, n = None):
#     """
#     This will plot the canonical result of the diffpop method for either whel or DMSO
#     :param df:
#     :return:
#     """
#     if gene_name not in dataframe['Genes'].values:
#         raise ValueError(f"Gene {gene_name} not found in the dataframe.")
#
#     treatment_columns = [col for col in dataframe.columns if treatment_type in col]
#     if not treatment_columns:
#         raise ValueError(f"No columns found for treatment type {treatment_type}.")
#
#     gene_row = dataframe[dataframe['Genes'] == gene_name]
#     x_labels = [col.split('-')[-1].replace('F', '') for col in treatment_columns]
#     y_values = gene_row[treatment_columns].values.flatten()
#
#     plt.figure(figsize=(10, 6))
#     plt.plot(x_labels, y_values, marker='o')
#     for i, txt in enumerate(y_values):
#         plt.annotate(txt, (x_labels[i], y_values[i]), textcoords="offset points", xytext=(0, 10), ha='center')
#     plt.xlabel('Fraction number')
#     plt.ylabel('Intensity')
#     plt.title(f'Protein Level for {gene_name} under {treatment_type} treatment run {n} ')
#     plt.grid(False)
#     plt.show()


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


def average_graph(dmso_df, whel_df, protein):
    """
    Plot the Average intensities with desired dfs.
    :param dmso_df:
    :param whel_df:
    :param protein:
    :return:
    """
    dmso = dmso_df[dmso_df["Genes"] == protein]
    whel = whel_df[whel_df["Genes"] == protein]

    dmso_values = dmso.drop(columns=['Genes']).values.flatten()
    whel_values = whel.drop(columns=['Genes']).values.flatten()

    # Get the x-axis labels (F1, F2, F3, ...)
    x_labels = ["2.62", "3.93", "5.91", "8.89", "13.38", "20.15", "30.38", "45.70", "100"]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(x_labels, dmso_values, marker='o', label='DMSO Dataset')
    plt.plot(x_labels, whel_values, marker='o', label='WHEL Dataset')

    plt.xlabel('Methanol Percentages')
    plt.ylabel(f'Log Transformed  Intensity of {protein} ')
    plt.title(f'Log Transformed Intensity for {protein} ')

    # plt.ylim(10,18)

    plt.legend()
    plt.grid(False)
    plt.savefig(f'img/Average of Log Transformed Intensity for {protein} .png')
    plt.close()


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


def average_dataframe(df):
    """
    Average the run of the dataframe for each fraction.
    :param df:
    :return:
    """
    genes = df['Genes']

    suffixes = df.columns.str.extract(r'-(F\d)')[0]

    averaged_columns = pd.DataFrame(index=df.index)

    for suffix in suffixes.dropna().unique():
        columns_to_average = df.columns[df.columns.str.contains(f'-{suffix}$')]
        prefix = df.columns[df.columns.str.contains(f'-{suffix}$')].str.extract(r'(\w+-)')[0][0]
        new_column_name = f'{prefix}{suffix}'
        averaged_columns[new_column_name] = df[columns_to_average].mean(axis=1)

    df_transformed = pd.concat([genes, averaged_columns], axis=1)

    return df_transformed


def log_df(df, string_col_name='Genes'):
    """
    Do a log transformation of the dataframe since it has a genes column.
    """
    string_column = df[string_col_name]
    numerical_columns = df.drop(columns=[string_col_name])

    adjusted_numerical_columns = numerical_columns
    log_transformed = np.log(adjusted_numerical_columns)

    log_transformed_df = pd.concat([string_column, log_transformed], axis=1)

    return log_transformed_df


def ks_test_between_runs(log_dmso, log_whel):
    """
    Test the KS score between every runs for DMSO and WHel in a combination. Ex: DMSO run1 vs Whel run1, whelrun2, whelrun3.
    :param log_dmso:
    :param log_whel:
    :return:
    """
    genes = log_dmso.index
    dmso_runs = sorted(set(col.split('-')[0] + '-' + col.split('-')[1] for col in log_dmso.columns))
    whel_runs = sorted(set(col.split('-')[0] + '-' + col.split('-')[1] for col in log_whel.columns))

    results = {'Genes': genes}

    for whel_run in whel_runs:
        for dmso_run in dmso_runs:
            column_name = f"{whel_run}_{dmso_run}"
            ks_results = []
            for gene in genes:
                whel_data = log_whel.loc[gene, [col for col in log_whel.columns if whel_run in col]].values
                dmso_data = log_dmso.loc[gene, [col for col in log_dmso.columns if dmso_run in col]].values
                ks_result = ks_2samp(whel_data, dmso_data)
                ks_results.append([ks_result.statistic, ks_result.pvalue])
            results[column_name] = ks_results

    return pd.DataFrame(results)


def ks_test_total(log_dmso, log_whel):
    """
    THe total score of the ks-test into a dataframe between two runs for all genes. THis is for the average scores.
    :param log_dmso:
    :param log_whel:
    :return:
    """

    genes = log_dmso.index
    ks_results = {'Genes': genes}
    results = []

    for gene in genes:
        dmso_data = log_dmso.loc[gene].values
        whel_data = log_whel.loc[gene].values
        ks_result = ks_2samp(dmso_data, whel_data)
        results.append([ks_result.statistic, ks_result.pvalue])

    ks_results['KS_Result'] = results
    return pd.DataFrame(ks_results)


def plot_ks_result_histogram(df, ks_column='KS_Result'):
    """
    Plot the histogram of the ks-score for downstream analysis.
    """
    df['KS_Result_First'] = df[ks_column].apply(lambda x: x[0])

    df_sorted = df.sort_values(by='KS_Result_First', ascending=False)

    unique_values = df['KS_Result_First'].unique()
    num_bins = len(unique_values)

    plt.figure(figsize=(10, 6))
    plt.hist(df['KS_Result_First'], bins=num_bins, edgecolor='black', color='orange')
    plt.title('Histogram of the KS Score')
    plt.xlabel('KS Score')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.show()

    return df_sorted


def plot_ks(df1, df2, gene):
    """
    Plot the KS plot of both DMSO and Wheldone on a specific gene.
    :param df1:
    :param df2:
    :param gene:
    :return:
    """
    data1 = df1.loc[gene].values
    data2 = df2.loc[gene].values

    ecdf1_x = np.sort(data1)
    ecdf2_x = np.sort(data2)
    ecdf1_y = np.arange(1, len(ecdf1_x) + 1) / len(ecdf1_x)
    ecdf2_y = np.arange(1, len(ecdf2_x) + 1) / len(ecdf2_x)

    ks_statistic, p_value = ks_2samp(data1, data2)

    plt.step(ecdf1_x, ecdf1_y, label=f'ECDF {gene} - DMSO', where='post')
    plt.step(ecdf2_x, ecdf2_y, label=f'ECDF {gene} - WHEL', where='post')

    max_diff_index = np.argmax(np.abs(ecdf1_y - ecdf2_y))
    plt.vlines(ecdf1_x[max_diff_index], ecdf2_y[max_diff_index], ecdf1_y[max_diff_index], colors='r', linestyle='dotted', label=f'KS Statistic: {ks_statistic:.4f}, p-value: {p_value:.4f}')

    plt.title(f'Kolmogorov-Smirnov Plot for {gene}')
    plt.xlabel('Log Transformed Intensity')
    plt.ylabel('Quantile')
    plt.legend()

    plt.savefig(f"KS_plot of {gene}.jpg")
    plt.close()


def plot_hist(df, gene):
    """
    plot th ehistogram of the data
    :param df:
    :param gene:
    :return:
    """
    data1 = df.loc[gene].values

    plt.hist(data1)

    plt.title(f'Histogram for {gene}')
    plt.xlabel('Log Transformed Intensity')
    plt.ylabel('Numbers')
    plt.legend()

    plt.savefig(f"Histogram for whel {gene}.jpg")
    plt.close()




def plot_gene_intensity(DMSO_df, whel_df, gene_name, normalize_mehod):
    """
    This is to plot the gene intensity with both the DMSO and Wheldone dataframe on a specific gene.
    :param DMSO_df:
    :param whel_df:
    :param gene_name:
    :return:
    """
    # Filter the dataframes to only include the specified gene
    DMSO_gene_data = DMSO_df[DMSO_df['Genes'] == gene_name]
    whel_gene_data = whel_df[whel_df['Genes'] == gene_name]

    # Reshape the DataFrames to long format
    DMSO_long = pd.melt(DMSO_gene_data, id_vars=['Genes'], var_name='Fraction', value_name='Intensity')
    whel_long = pd.melt(whel_gene_data, id_vars=['Genes'], var_name='Fraction', value_name='Intensity')

    # Extract run and fraction information
    DMSO_long[['Treatment', 'Run', 'Fraction']] = DMSO_long['Fraction'].str.extract(r'(DMSO)-n(\d+)-F(\d+)')
    whel_long[['Treatment', 'Run', 'Fraction']] = whel_long['Fraction'].str.extract(r'(whel)-n(\d+)-F(\d+)')

    # Combine the two DataFrames
    combined_df = pd.concat([DMSO_long, whel_long])

    # Convert relevant columns to numeric
    combined_df['Run'] = pd.to_numeric(combined_df['Run'])
    combined_df['Fraction'] = pd.to_numeric(combined_df['Fraction'])

    plt.figure(figsize=(12, 8))
    sns.lineplot(data=combined_df, x='Fraction', y='Intensity', hue='Treatment', style='Treatment', markers=True,
                 errorbar ='sd', err_style='band')

    x_labels = [f"F{i}" for i in range(0,9)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)

    # plt.ylim(-0.2, 1)

    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_mehod} normalized Intensity')
    plt.title(f'{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run')
    plt.legend(title='Treatment')
    plt.grid(True)

    # Save the plot
    plt.savefig(f"img/{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run.png")
    plt.close()


def one_to_five_and_to_nine_average(df):
    """
    Check the difference between the Fraction 1-5 and Fraction 5-9. We only want the runs with higher average in the later
    fraction.
    :param df:
    :return:
    """
    df_copy = df.copy()

    df_copy['avg_1_to_5'] = df.iloc[:, 1:5].mean(axis=1)

    df_copy['avg_5_to_9'] = df.iloc[:, 5:9].mean(axis=1)

    df_copy['difference'] = df_copy['avg_5_to_9'] - df_copy['avg_1_to_5']

    df_result = df_copy.iloc[:, [0, -3, -2, -1]]

    df_filtered = df_result[df_result['avg_5_to_9'] > df_result['avg_1_to_5']]

    df_sorted = df_filtered.sort_values(by='difference', ascending=False)

    return df_sorted


def z_score_normalization(column):
    return (column - column.mean()) / column.std()


def z_norm_df(df):
    """
    Nramalize the data with z-normalization method
    :param df:
    :return:
    """
    df_genes = df["Genes"]
    numerical_columns = df.drop(columns=["Genes"])
    df_standardized = numerical_columns.apply(z_score_normalization)
    df_standardized.insert(0, 'Genes', df_genes)

    return df_standardized


def nsaf_normalization(df):
    """
    Perform NSAF normalization to the given df for all columns despite the first ("Genes") and the last ("Length") column
    :param df:
    :return:
    """
    first_column = df.columns[0]
    run_columns = df.columns[1:-1]
    length_column = df.columns[-1]
    nsaf_df = pd.DataFrame(df[first_column])

    # NSAF for each column.
    for run in run_columns:
        df[f'SAF_{run}'] = df[run] / df[length_column]
        total_saf = df[f'SAF_{run}'].sum()
        nsaf_df[f'{run}'] = df[f'SAF_{run}'] / total_saf

    return nsaf_df


def tic_normalization(df):
    """
    Perform TIC normalization on the given DataFrame.
    """
    # Extract the columns to normalize (the first column is 'Genes')
    intensity_columns = df.columns[1:]

    # Calculate the total ion current for each sample
    tic = df[intensity_columns].sum(axis=0)

    # Perform TIC normalization
    df_normalized = df.copy()
    df_normalized[intensity_columns] = df[intensity_columns].div(tic, axis=1)

    return df_normalized


def median_normalization(df):
    """
    Perform median normalization on the given DataFrame.
    """
    # Extract the columns to normalize (the first column is 'Genes')
    intensity_columns = df.columns[1:]

    # Calculate the median intensity for each sample
    median_intensities = df[intensity_columns].median(axis=0)

    # Perform median normalization
    df_normalized = df.copy()
    df_normalized[intensity_columns] = df[intensity_columns].div(median_intensities, axis=1)

    return df_normalized


def quantile_normalization(df):
    """
    Perform quantile normalization on the given DataFrame.
    """
    # Extract the columns to normalize (the first column is 'Genes')
    genes_column = df.iloc[:, 0]
    intensity_columns = df.columns[1:]

    # Initialize the QuantileTransformer
    transformer = QuantileTransformer(output_distribution='uniform', random_state=0)

    # Perform quantile normalization
    df_normalized = df.copy()
    df_normalized[intensity_columns] = transformer.fit_transform(df[intensity_columns])

    df_normalized.iloc[:, 0] = genes_column

    return df_normalized


ccp = pd.read_excel("ccp.xlsx", "Atlas")
protein_length = pd.read_csv("protein.tsv", delimiter= "\t")
protein_length_2 = pd.read_csv("protein_rerun.tsv", delimiter="\t")
protein_length_df = protein_length[["Genes","Length"]]
protein_length_df_2 = protein_length_2[["Genes","Length"]]
union_protein_length = pd.concat([protein_length_df, protein_length_df_2])

df_new = pd.read_csv("new_dataset.csv")
df_dmso, df_whel = separate_dataframe(df_new)
filled_dmso = fill_na_with_half_min(df_dmso).dropna()
filled_whel = fill_na_with_half_min(df_whel).dropna()

# Get all the genes for each run and Union them
dmso_gene = set(filled_dmso["Genes"])
whel_gene = set(filled_whel["Genes"])
inter_genes = list(dmso_gene & whel_gene)

# Subset shared genes
dmso_shared = filled_dmso[filled_dmso["Genes"].isin(inter_genes)]
whel_shared = filled_whel[filled_whel["Genes"].isin(inter_genes)]

if __name__ == '__main__':

    print("Hello World!")

    # Read the files as dataframe
    # ccp = pd.read_excel("ccp.xlsx", "Atlas")
    # protein_length = pd.read_csv("protein.tsv", delimiter= "\t")
    # protein_length_2 = pd.read_csv("protein_rerun.tsv", delimiter="\t")
    # protein_length_df = protein_length[["Genes","Length"]]
    # protein_length_df_2 = protein_length_2[["Genes","Length"]]
    # union_protein_length = pd.concat([protein_length_df, protein_length_df_2])
    #
    # df_new = pd.read_csv("new_dataset.csv")
    # df_dmso, df_whel = separate_dataframe(df_new)
    # filled_dmso = fill_na_with_half_min(df_dmso).dropna()
    # filled_whel = fill_na_with_half_min(df_whel).dropna()
    #
    # # Get all the genes for each run and Union them
    # dmso_gene = set(filled_dmso["Genes"])
    # whel_gene = set(filled_whel["Genes"])
    # inter_genes = list(dmso_gene & whel_gene)
    #
    # # Subset shared genes
    # dmso_shared = filled_dmso[filled_dmso["Genes"].isin(inter_genes)]
    # whel_shared = filled_whel[filled_whel["Genes"].isin(inter_genes)]


    """ 
    # The NSFA normalization
    """
    # Add the protein length for NSFA method
    nsfa_df_dmso = pd.merge(dmso_shared, union_protein_length, on = "Genes", how = "left").drop_duplicates(subset=['Genes'])
    nsfa_df_whel = pd.merge(whel_shared, union_protein_length, on = "Genes", how = "left").drop_duplicates(subset=['Genes'])

    nsaf_normalized_dmso = nsaf_normalization(nsfa_df_dmso)
    nsaf_normalized_whel = nsaf_normalization(nsfa_df_whel)
    


    """
    The TIC normalization
    """
    tic_normalized_dmso = tic_normalization(dmso_shared)
    tic_normalized_whel = tic_normalization(whel_shared)

    """
    Median normalization
    """
    median_normalized_dmso = median_normalization(dmso_shared)
    median_normalized_whel = median_normalization(whel_shared)

    """
    Quantile normalization
    """
    quantile_normalized_dmso = quantile_normalization(dmso_shared)
    quantile_normalized_whel = quantile_normalization(whel_shared)

    # Plot nsaf normalizations
    for i in union_set:
        plot_gene_intensity(nsaf_normalized_dmso, nsaf_normalized_whel, i, "NSAF")

    """
    This part is the Z-Norm normalization
    
    # # Z-Norm Normalization
    # dmso_norm = z_norm_df(dmso_shared)
    # whel_norm = z_norm_df(whel_shared)
    # 
    # # Average the normalized runs
    # avg_norm_dmso = average_dataframe(dmso_norm)
    # avg_norm_whel = average_dataframe(whel_norm)
    # 
    # plot_gene_intensity(dmso_norm, whel_norm, "AURKA")
    # plot_gene_intensity(dmso_norm, whel_norm, "AURKB")
    # plot_gene_intensity(dmso_norm, whel_norm, "PRC1")
    # plot_gene_intensity(dmso_norm, whel_norm, "KIF11")
    # plot_gene_intensity(dmso_norm, whel_norm, "CCNB1")
    # plot_gene_intensity(dmso_norm, whel_norm, "TACC3")
    # 
    # average_graph(avg_norm_dmso, avg_norm_whel, "AURKA")
    # average_graph(avg_norm_dmso, avg_norm_whel, "AURKB")
    # average_graph(avg_norm_dmso, avg_norm_whel, "PRC1")
    # average_graph(avg_norm_dmso, avg_norm_whel, "KIF11")
    # average_graph(avg_norm_dmso, avg_norm_whel, "CCNB1")
    # average_graph(avg_norm_dmso, avg_norm_whel, "TACC3")
    """


    #

    # Log Transformed
    # log_dmso = log_df(dmso_shared)
    # log_whel = log_df(whel_shared)
    #
    # subset_dmso = log_dmso.set_index("Genes")
    #
    # subset_whel = log_whel.set_index("Genes")
    #
    # # total_res = ks_test_total(subset_dmso, subset_whel)
    #
    # # sorted_df = plot_ks_result_histogram(total_res)
    #
    # # join_df = pd.merge(sorted_df, ccp, on ="Genes", how = "inner")
    #
    # # result_df = ks_test_between_runs(subset_dmso, subset_whel)
    #
    # avg_dmso = average_dataframe(dmso_shared)
    # avg_whel = average_dataframe(whel_shared)
    #
    # log_average_dmso = log_df(avg_dmso)
    # log_average_whel = log_df(avg_whel)
    #
    # avg_log_dmso = average_dataframe(log_dmso)
    # avg_log_whel = average_dataframe(log_whel)



    #
    # filled_whel.to_excel("filled_whel.xlsx", index = False)
    # filled_dmso.to_excel("filled_dmso.xlsx", index = False)
    #
    # # Plot the average of run1,2,3 of the protein for DiffPoP
    # average_graph(log_average_dmso, log_average_whel, "AURKA")
    # average_graph(log_average_dmso, log_average_whel, "AURKB")
    # average_graph(log_average_dmso, log_average_whel, "PRC1")
    # average_graph(log_average_dmso, log_average_whel, "KIF11")
    # average_graph(log_average_dmso, log_average_whel, "CCNB1")
    # average_graph(log_average_dmso, log_average_whel, "TACC3")
    #
    # plot_ks(subset_dmso,subset_whel, "AURKA")
    # plot_ks(subset_dmso,subset_whel, "AURKB")
    # plot_ks(subset_dmso,subset_whel, "PRC1")
    # plot_ks(subset_dmso,subset_whel, "KIF11")
    # plot_ks(subset_dmso,subset_whel, "CCNB1")
    # plot_ks(subset_dmso,subset_whel, "TACC3")

    # plot_hist(subset_whel, "CCNB1")

    # plot_gene_intensity(log_dmso, log_whel, "AURKA")
    # plot_gene_intensity(log_dmso, log_whel, "AURKB")
    # plot_gene_intensity(log_dmso, log_whel, "PRC1")
    # plot_gene_intensity(log_dmso, log_whel, "KIF11")
    # plot_gene_intensity(log_dmso, log_whel, "CCNB1")
    # plot_gene_intensity(log_dmso, log_whel, "TACC3")
    #
    # average_graph(avg_log_dmso, avg_log_whel, "AURKA")
    # average_graph(avg_log_dmso, avg_log_whel, "AURKB")
    # average_graph(avg_log_dmso, avg_log_whel, "PRC1")
    # average_graph(avg_log_dmso, avg_log_whel, "KIF11")
    # average_graph(avg_log_dmso, avg_log_whel, "CCNB1")
    # average_graph(avg_log_dmso, avg_log_whel, "TACC3")
    #
    # a_of_a_dmso = one_to_five_and_to_nine_average(avg_log_dmso)
    # a_of_a_whel = one_to_five_and_to_nine_average(avg_log_whel)

    #
    # merged_aoa = pd.merge(a_of_a_dmso, a_of_a_whel, on= "Genes", how= "inner")

    # for i in merged_aoa["Genes"]:
    #     plot_gene_intensity(log_dmso, log_whel, i)
    #     average_graph(avg_log_dmso, avg_log_whel, i)

    # KS_df = pd.read_excel("Sorted_KS-SCORE_with_threshold_0.1_DMSO_vs_whel.xlsx")
    # AURKA_df = KS_df[KS_df["Genes"].isin(AURKA)]
    # AURKB_df = KS_df[KS_df["Genes"].isin(AURKB)]
    # KIF11_df = KS_df[KS_df["Genes"].isin(KIF11)]
    # PRC1_df = KS_df[KS_df["Genes"].isin(PRC1)]
    # CCNB1_df = KS_df[KS_df["Genes"].isin(CCNB1)]
    #
    # AURKA_df.to_excel("AURKA.xlsx", index = False)
    # AURKB_df.to_excel("AURKB.xlsx", index = False)
    # KIF11_df.to_excel("KIF11.xlsx", index = False)
    # PRC1_df.to_excel("PRC1.xlsx", index = False)
    # CCNB1_df.to_excel("CCNB1.xlsx", index = False)
    #
    # AURKA_list = AURKA_df["Genes"].tolist()
    # AURKB_list = AURKB_df["Genes"].tolist()
    # KIF11_list = KIF11_df["Genes"].tolist()
    # PRC1_list = PRC1_df["Genes"].tolist()
    # CCNB1_list = CCNB1_df["Genes"].tolist()
    #
    # union_set = list(union_lists_to_set(AURKA_list, AURKB_list, KIF11_list, PRC1_list, CCNB1_list))
    #
    # for i in union_set:
    #     plot_gene_intensity(log_dmso, log_whel, i)
    #     average_graph(avg_log_dmso, avg_log_whel, i)
    #

    # #This is for the positive runs
    # averaged_run = pd.read_excel("Positive_fraction.xlsx")
    # dmso_df, whel_df = separate_dataframe(averaged_run)
    #
    # average_graph(dmso_df, whel_df, "MYCBP")
    # average_graph(dmso_df, whel_df, "MYCBP2")
    # average_graph(dmso_df, whel_df, "AURKA")
    # average_graph(dmso_df, whel_df, "AURKB")
    # average_graph(dmso_df, whel_df, "TPX2")


    # #This is for the averaged runs
    # averaged_run = pd.read_csv("Averaegd_runs.csv")
    # dmso_df, whel_df = separate_dataframe(averaged_run)
    #
    # average_graph(dmso_df, whel_df, "MYCBP")
    # average_graph(dmso_df, whel_df, "MYCBP2")
    # average_graph(dmso_df, whel_df, "AURKA")
    # average_graph(dmso_df, whel_df, "AURKB")
    # average_graph(dmso_df, whel_df, "TPX2")



    # # This is for individual runs
    # experiment = pd.read_excel("With_R1_3_Fractioned.xlsx")
    # columns_selection_df = experiment.iloc[:1]
    # gene_whel_runs, gene_dmso_runs = whel_dmso_subset(columns_selection_df)
    # g1 = gene_dmso_runs["1"]
    # g2 = gene_dmso_runs["2"]
    # g3 = gene_dmso_runs["3"]
    # w1 = gene_whel_runs["1"]
    # w2 = gene_whel_runs["2"]
    # w3 = gene_whel_runs["3"]

