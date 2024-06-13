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


def nsaf_method(df):
    exclude_columns = ['Genes', 'Protein Length']
    columns_to_divide = [col for col in df.columns if col not in exclude_columns]

    # Perform the division
    for col in columns_to_divide:
        df[col] = df[col] / df['Protein Length']

    sum_safs = df[columns_to_divide].sum()

    # Calculate NSAF for each column
    for col in columns_to_divide:
        df[col] = df[col] / sum_safs[col]

    return df


# def nsaf_normalization(df_f):
#     """
#     Perform NSAF normalization to the given df for all columns despite the first ("Genes") and the last ("Length") column
#     :param df:
#     :return:
#     """
#     df = df_f.copy()
#     first_column = df.columns[0]
#     run_columns = df.columns[1:-2]
#     length_column = 'Length'
#     total_spec_column = "Total Spectral Count"
#
#     # Create a new DataFrame to store NSAF results
#     nsaf_df = pd.DataFrame(df[first_column])
#
#     # Calculate SAF for each protein
#     df['SAF'] = df[total_spec_column] / df[length_column]
#
#     # Calculate the total SAF for the sample
#     total_saf = df['SAF'].sum()
#
#     # Calculate NSAF for each protein
#     df['NSAF'] = df['SAF'] / total_saf
#
#     # Multiply NSAF by the intensity for each run
#     for run in run_columns:
#         nsaf_intensity_column = f'{run}'
#         df[nsaf_intensity_column] = df['NSAF'] * df[run]
#         nsaf_df[nsaf_intensity_column] = df[nsaf_intensity_column]
#
#     return nsaf_df


def nsaf_func(df_f):
    df = df_f.copy()
    first_column = df.columns[0]
    run_columns = df.columns[1:-2]
    length_column = 'Length'
    total_spec_column = "Razor Spectral Count"

    # Create a new DataFrame to store NSAF results
    nsaf_df = pd.DataFrame(df[first_column])

    # Calculate SAF for each protein
    df['SAF'] = df[total_spec_column] / df[length_column]

    # Calculate the total SAF for the sample
    total_saf = df['SAF'].sum()

    # Calculate NSAF for each protein
    df['NSAF'] = df['SAF'] / total_saf

    nsaf_df["NSAF"] = df["NSAF"]
    nsaf_df["Razor Spectral Count"] = df[total_spec_column]


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


if __name__ == '__main__':
    print("Hello")