from read_my_file import *


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

    return gene_whel_data, gene_dmso_data

def graphing_methods(df):
    """
    This will plot the canonical result of the diffpop method for either whel or DMSO
    :param df:
    :return:
    """

if __name__ == '__main__':

    experiment = pd.read_excel("With_R1.xlsx")
    columns_selection_df = experiment.iloc[:, [3] + list(range(5, experiment.shape[1]))]
    gene_whel_df, gene_dmso_df = whel_dmso_subset(columns_selection_df)