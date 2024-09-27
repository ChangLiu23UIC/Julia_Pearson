import seaborn as sns
from sklearn.preprocessing import QuantileTransformer
import pandas as pd
import numpy as np


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


def var_stab_normalization(df):
    """
    Perform Variance Stabilization Normalization (VSN) on the given DataFrame.
    """
    # Extract the columns to normalize (the first column is 'Genes')
    genes_column = df.iloc[:, 0]
    intensity_columns = df.columns[1:]

    # Apply log transformation to stabilize variance
    df_log_transformed = np.log2(df[intensity_columns] + 1)

    # Scale data to have mean 0 and standard deviation 1
    df_mean = df_log_transformed.mean()
    df_std = df_log_transformed.std()
    df_vsn_normalized = (df_log_transformed - df_mean) / df_std

    # Reattach the genes column
    df_vsn_normalized.insert(0, df.columns[0], genes_column)

    return df_vsn_normalized