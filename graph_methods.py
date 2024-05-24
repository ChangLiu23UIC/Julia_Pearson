from read_my_file import *
import matplotlib.pyplot as plt


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

def graphing_methods(dataframe, gene_name):
    """
    This will plot the canonical result of the diffpop method for either whel or DMSO
    :param df:
    :return:
    """
    if gene_name not in dataframe['Genes'].values:
        raise ValueError(f"Gene {gene_name} not found in the dataframe.")

    gene_row = dataframe[dataframe['Genes'] == gene_name]
    x_labels = [col.split('-')[-1].replace('F', '') for col in dataframe.columns if 'DMSO' in col]
    y_values = gene_row[[col for col in dataframe.columns if 'DMSO' in col]].values.flatten()

    plt.figure(figsize=(10, 6))
    plt.plot(x_labels, y_values, marker='o')
    for i, txt in enumerate(y_values):
        plt.annotate(txt, (x_labels[i], y_values[i]), textcoords="offset points", xytext=(0, 10), ha='center')
    plt.xlabel('F number')
    plt.ylabel('Intensity')
    plt.title(f'Gene Expression for {gene_name}')
    plt.grid(True)
    plt.show()


if __name__ == '__main__':

    experiment = pd.read_excel("With_R1.xlsx")
    columns_selection_df = experiment.iloc[:, [3] + list(range(5, experiment.shape[1]))]
    columns_selection_df = columns_selection_df.fillna(0)
    gene_whel_runs, gene_dmso_runs = whel_dmso_subset(columns_selection_df)
    g1 = gene_dmso_runs["1"]
    g2 = gene_dmso_runs["2"]
    g3 = gene_dmso_runs["3"]
    w1 = gene_whel_runs["1"]
    w2 = gene_whel_runs["2"]
    w3 = gene_whel_runs["3"]