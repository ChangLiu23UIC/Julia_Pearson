import seaborn as sns
from sklearn.preprocessing import QuantileTransformer
import pandas as pd


def z_normalization(df):
    """
    Perform Z-score normalization on the given DataFrame.
    """
    # Extract the columns to normalize (the first column is 'Genes')
    genes_column = df.iloc[:, 0]
    intensity_columns = df.columns[1:]

    # Perform Z-score normalization
    df_normalized = df.copy()
    df_normalized[intensity_columns] = (df[intensity_columns] - df[intensity_columns].mean()) / df[
        intensity_columns].std()

    df_normalized.iloc[:, 0] = genes_column

    return df_normalized


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

if __name__ == '__main__':
    print("Hello")