from read_my_file import *
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import numpy as np


def whel_dmso_subset(df):
    """
    Subset the dataframe into gene + whel-xxx and gene + dmso-xxx
    :param df:
    :return:
    """
    gene_column = [col for col in df.columns if col == 'Genes']
    whel_columns = [col for col in df.columns if col.startswith('whel-')]
    dmso_columns = [col for col in df.columns if col.startswith('DMSO-')]

    gene_whel_data = df[gene_column + whel_columns]
    gene_dmso_data = df[gene_column + dmso_columns]

    dmso_runs = subset_by_run(gene_dmso_data)
    whel_runs = subset_by_run(gene_whel_data)

    return whel_runs, dmso_runs


def subset_by_run(df):
    """

    :param df:
    :return:
    """
    genes_column = df[['Genes']]

    x_values = set(col.split('-')[1][1:] for col in df.columns if '-' in col and col != 'Genes')

    subsets = {}

    for x in x_values:
        columns_with_x = [col for col in df.columns if f'-r{x}-' in col or f'-n{x}-' in col]
        columns_with_x_sorted = sorted(columns_with_x, key=lambda col: int(col.split('-F')[-1]))
        subset = pd.concat([genes_column, df[columns_with_x_sorted]], axis=1)
        subsets[x] = subset
    return subsets

def graphing_methods(dataframe, gene_name,treatment_type, n = None):
    """
    This will plot the canonical result of the diffpop method for either whel or DMSO
    :param df:
    :return:
    """
    if gene_name not in dataframe['Genes'].values:
        raise ValueError(f"Gene {gene_name} not found in the dataframe.")

    treatment_columns = [col for col in dataframe.columns if treatment_type in col]
    if not treatment_columns:
        raise ValueError(f"No columns found for treatment type {treatment_type}.")

    gene_row = dataframe[dataframe['Genes'] == gene_name]
    x_labels = [col.split('-')[-1].replace('F', '') for col in treatment_columns]
    y_values = gene_row[treatment_columns].values.flatten()

    plt.figure(figsize=(10, 6))
    plt.plot(x_labels, y_values, marker='o')
    for i, txt in enumerate(y_values):
        plt.annotate(txt, (x_labels[i], y_values[i]), textcoords="offset points", xytext=(0, 10), ha='center')
    plt.xlabel('Fraction number')
    plt.ylabel('Intensity')
    plt.title(f'Protein Level for {gene_name} under {treatment_type} treatment run {n} ')
    plt.grid(False)
    plt.show()


def separate_dataframe(df):
    genes_col = df['Genes']

    dmso_cols = [col for col in df.columns if col.startswith('DMSO-')]
    whel_cols = [col for col in df.columns if col.startswith('whel-')]

    df_dmso = pd.DataFrame({'Genes': genes_col})
    df_dmso = pd.concat([df_dmso, df[dmso_cols]], axis=1)

    df_whel = pd.DataFrame({'Genes': genes_col})
    df_whel = pd.concat([df_whel, df[whel_cols]], axis=1)

    return df_dmso, df_whel


def average_graph(dmso_df, whel_df, protein):
    dmso = dmso_df[dmso_df["Genes"] == protein]
    whel = whel_df[whel_df["Genes"] == protein]

    dmso_values = dmso.drop(columns=['Genes']).values.flatten()
    whel_values = whel.drop(columns=['Genes']).values.flatten()

    # Get the x-axis labels (F1, F2, F3, ...)
    x_labels = [f'F{i}' for i in range(1, len(dmso_values) + 1)]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(x_labels, dmso_values, marker='o', label='DMSO Dataset')
    plt.plot(x_labels, whel_values, marker='o', label='WHEL Dataset')

    plt.xlabel('Fraction')
    plt.ylabel(f'Absolute Intensity of {protein} protein')
    plt.title(f'Absolute for {protein} protein')
    plt.legend()
    plt.grid(False)
    plt.savefig(f'Absolute Intensity for {protein} protein.png')
    plt.close()


def fill_na_with_half_min(df):
    """
    """
    categorical_col = df.iloc[:, 0]
    numeric_df = df.iloc[:, 1:]

    filled_numeric_df = numeric_df.apply(lambda row: row.fillna(row.min(skipna=True)/2), axis=1)
    df_filled = pd.concat([categorical_col, filled_numeric_df], axis=1)

    return df_filled


def transform_dataframe(df):
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
    """
    string_column = df[string_col_name]
    numerical_columns = df.drop(columns=[string_col_name])

    adjusted_numerical_columns = numerical_columns
    log_transformed = np.log(adjusted_numerical_columns)

    log_transformed_df = pd.concat([string_column, log_transformed], axis=1)

    return log_transformed_df


def ks_test_between_runs(log_dmso, log_whel):
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


if __name__ == '__main__':

    print("Hello World!")

    df_new = pd.read_csv("new_dataset.csv")
    df_dmso, df_whel = separate_dataframe(df_new)
    filled_dmso = fill_na_with_half_min(df_dmso).dropna()
    filled_whel = fill_na_with_half_min(df_whel).dropna()

    dmso_gene = set(filled_dmso["Genes"])
    whel_gene = set(filled_whel["Genes"])

    inter_genes = list(dmso_gene & whel_gene)

    dmso_shared = filled_dmso[filled_dmso["Genes"].isin(inter_genes)]
    whel_shared = filled_whel[filled_whel["Genes"].isin(inter_genes)]

    log_dmso = log_df(dmso_shared)
    log_whel = log_df(whel_shared)

    subset_dmso = log_dmso.set_index("Genes")

    subset_whel = log_whel.set_index("Genes")


    result_df = ks_test_between_runs(subset_dmso, subset_whel)




    # avg_dmso = transform_dataframe(filled_dmso)
    # avg_whel = transform_dataframe(filled_whel)
    #
    # filled_whel.to_excel("filled_whel.xlsx", index = False)
    # filled_dmso.to_excel("filled_dmso.xlsx", index = False)
    #
    #
    # average_graph(avg_dmso, avg_whel, "MYCBP")
    # average_graph(avg_dmso, avg_whel, "MYCBP2")
    # average_graph(avg_dmso, avg_whel, "AURKA")
    # average_graph(avg_dmso, avg_whel, "AURKB")
    # average_graph(avg_dmso, avg_whel, "TPX2")
























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

