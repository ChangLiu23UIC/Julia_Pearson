import seaborn as sns
from sklearn.preprocessing import QuantileTransformer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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


def all_norm_methods(df):
    df_z = z_normalization(df)
    df_tic = tic_normalization(df)
    df_med = median_normalization(df)
    df_quan = quantile_normalization(df)
    df_var = var_stab_normalization(df)

    return df_z, df_tic, df_med, df_quan, df_var


def plot_gene_intensity_DP(DMSO_df, peng_df, gene_name, normalize_mehod):
    """
    This is to plot the gene intensity with both the DMSO and Wheldone dataframe on a specific gene.
    :param DMSO_df:
    :param whel_df:
    :param gene_name:
    :return:
    """
    # Filter the dataframes to only include the specified gene
    DMSO_gene_data = DMSO_df[DMSO_df['Gene'] == gene_name]
    # cipr_gene_data = cipr_df[cipr_df['Gene'] == gene_name]
    # moxi_gene_data = moxi_df[moxi_df['Gene'] == gene_name]
    peng_gene_data = peng_df[peng_df['Gene'] == gene_name]


    # Reshape the DataFrames to long format
    DMSO_long = pd.melt(DMSO_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    # cipr_long = pd.melt(cipr_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    # moxi_long = pd.melt(moxi_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    peng_long = pd.melt(peng_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')


    # Extract run and fraction information
    DMSO_long[['Treatment', 'Fraction']] = DMSO_long['Fraction'].str.extract(r'(DMSO)-F(\d+)')
    # cipr_long[['Treatment', 'Fraction']] = cipr_long['Fraction'].str.extract(r'(Cipro)-F(\d+)')
    # moxi_long[['Treatment', 'Fraction']] = moxi_long['Fraction'].str.extract(r'(Moxi)-F(\d+)')
    peng_long[['Treatment', 'Fraction']] = peng_long['Fraction'].str.extract(r'(PenG)-F(\d+)')


    # Combine the two DataFrames
    combined_df = pd.concat([DMSO_long, peng_long])

    # Convert relevant columns to numeric
    combined_df['Fraction'] = pd.to_numeric(combined_df['Fraction'])

    plt.figure(figsize=(12, 8))
    sns.lineplot(data=combined_df, x='Fraction', y='Intensity', hue='Treatment', style='Treatment', markers=True,
                 errorbar ='sd', err_style='band')

    x_labels = [f"F{i}" for i in range(1,10)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)

    # plt.ylim(-0.2, 1)

    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_mehod} normalized Intensity')
    plt.title(f'{normalize_mehod} normalized Intensity for {gene_name} for each')
    plt.legend(title='Treatment')
    plt.grid(True)
    plt.show()

    # Save the plot
    plt.savefig(f"Result/{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run.png")


def plot_gene_intensity_CM(DMSO_df, cipr_df, gene_name, normalize_mehod):
    """
    This is to plot the gene intensity with both the DMSO and Wheldone dataframe on a specific gene.
    :param DMSO_df:
    :param whel_df:
    :param gene_name:
    :return:
    """
    # Filter the dataframes to only include the specified gene
    DMSO_gene_data = DMSO_df[DMSO_df['Gene'] == gene_name]
    cipr_gene_data = cipr_df[cipr_df['Gene'] == gene_name]
    # moxi_gene_data = moxi_df[moxi_df['Gene'] == gene_name]
    # peng_gene_data = peng_df[peng_df['Gene'] == gene_name]


    # Reshape the DataFrames to long format
    DMSO_long = pd.melt(DMSO_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    cipr_long = pd.melt(cipr_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    # moxi_long = pd.melt(moxi_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    # peng_long = pd.melt(peng_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')


    # Extract run and fraction information
    DMSO_long[['Treatment', 'Fraction']] = DMSO_long['Fraction'].str.extract(r'(DMSO)-F(\d+)')
    cipr_long[['Treatment', 'Fraction']] = cipr_long['Fraction'].str.extract(r'(Cipro)-F(\d+)')
    # moxi_long[['Treatment', 'Fraction']] = moxi_long['Fraction'].str.extract(r'(Moxi)-F(\d+)')
    # peng_long[['Treatment', 'Fraction']] = peng_long['Fraction'].str.extract(r'(PenG)-F(\d+)')


    # Combine the two DataFrames
    combined_df = pd.concat([DMSO_long, cipr_long])

    # Convert relevant columns to numeric
    combined_df['Fraction'] = pd.to_numeric(combined_df['Fraction'])

    plt.figure(figsize=(12, 8))
    sns.lineplot(data=combined_df, x='Fraction', y='Intensity', hue='Treatment', style='Treatment', markers=True,
                 errorbar ='sd', err_style='band')

    x_labels = [f"F{i}" for i in range(1,10)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)

    # plt.ylim(-0.2, 1)

    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_mehod} normalized Intensity')
    plt.title(f'{normalize_mehod} normalized Intensity for {gene_name} for each')
    plt.legend(title='Treatment')
    plt.grid(True)
    plt.show()

    # Save the plot
    plt.savefig(f"Result/{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run.png")

def plot_gene_intensity_DM(DMSO_df, moxi_df, gene_name, normalize_mehod):
    """
    This is to plot the gene intensity with both the DMSO and Wheldone dataframe on a specific gene.
    :param DMSO_df:
    :param whel_df:
    :param gene_name:
    :return:
    """
    # Filter the dataframes to only include the specified gene
    DMSO_gene_data = DMSO_df[DMSO_df['Gene'] == gene_name]
    # cipr_gene_data = cipr_df[cipr_df['Gene'] == gene_name]
    moxi_gene_data = moxi_df[moxi_df['Gene'] == gene_name]
    # peng_gene_data = peng_df[peng_df['Gene'] == gene_name]


    # Reshape the DataFrames to long format
    DMSO_long = pd.melt(DMSO_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    # cipr_long = pd.melt(cipr_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    moxi_long = pd.melt(moxi_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    # peng_long = pd.melt(peng_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')


    # Extract run and fraction information
    DMSO_long[['Treatment', 'Fraction']] = DMSO_long['Fraction'].str.extract(r'(DMSO)-F(\d+)')
    # cipr_long[['Treatment', 'Fraction']] = cipr_long['Fraction'].str.extract(r'(Cipro)-F(\d+)')
    moxi_long[['Treatment', 'Fraction']] = moxi_long['Fraction'].str.extract(r'(Moxi)-F(\d+)')
    # peng_long[['Treatment', 'Fraction']] = peng_long['Fraction'].str.extract(r'(PenG)-F(\d+)')


    # Combine the two DataFrames
    combined_df = pd.concat([DMSO_long, moxi_long])

    # Convert relevant columns to numeric
    combined_df['Fraction'] = pd.to_numeric(combined_df['Fraction'])

    plt.figure(figsize=(12, 8))
    sns.lineplot(data=combined_df, x='Fraction', y='Intensity', hue='Treatment', style='Treatment', markers=True,
                 errorbar ='sd', err_style='band')

    x_labels = [f"F{i}" for i in range(1,10)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)

    # plt.ylim(-0.2, 1)

    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_mehod} normalized Intensity')
    plt.title(f'{normalize_mehod} normalized Intensity for {gene_name} for each')
    plt.legend(title='Treatment')
    plt.grid(True)
    plt.show()

    # Save the plot
    plt.savefig(f"Result/{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run.png")


def plot_gene_intensity_all(DMSO_df,cipr_df, moxi_df, peng_df, gene_name, normalize_mehod):
    """
    This is to plot the gene intensity with both the DMSO and Wheldone dataframe on a specific gene.
    :param DMSO_df:
    :param whel_df:
    :param gene_name:
    :return:
    """
    # Filter the dataframes to only include the specified gene
    DMSO_gene_data = DMSO_df[DMSO_df['Gene'] == gene_name]
    cipr_gene_data = cipr_df[cipr_df['Gene'] == gene_name]
    moxi_gene_data = moxi_df[moxi_df['Gene'] == gene_name]
    peng_gene_data = peng_df[peng_df['Gene'] == gene_name]


    # Reshape the DataFrames to long format
    DMSO_long = pd.melt(DMSO_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    cipr_long = pd.melt(cipr_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    moxi_long = pd.melt(moxi_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')
    peng_long = pd.melt(peng_gene_data, id_vars=['Gene'], var_name='Fraction', value_name='Intensity')


    # Extract run and fraction information
    DMSO_long[['Treatment', 'Fraction']] = DMSO_long['Fraction'].str.extract(r'(DMSO)-F(\d+)')
    cipr_long[['Treatment', 'Fraction']] = cipr_long['Fraction'].str.extract(r'(Cipro)-F(\d+)')
    moxi_long[['Treatment', 'Fraction']] = moxi_long['Fraction'].str.extract(r'(Moxi)-F(\d+)')
    peng_long[['Treatment', 'Fraction']] = peng_long['Fraction'].str.extract(r'(PenG)-F(\d+)')


    # Combine the two DataFrames
    combined_df = pd.concat([DMSO_long, cipr_long, moxi_long, peng_long])

    # Convert relevant columns to numeric
    combined_df['Fraction'] = pd.to_numeric(combined_df['Fraction'])

    plt.figure(figsize=(12, 8))
    sns.lineplot(data=combined_df, x='Fraction', y='Intensity', hue='Treatment', style='Treatment', markers=True,
                 errorbar ='sd', err_style='band')

    x_labels = [f"F{i}" for i in range(1,10)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)

    # plt.ylim(-0.2, 1)

    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_mehod} normalized Intensity')
    plt.title(f'{normalize_mehod} normalized Intensity for {gene_name} for each')
    plt.legend(title='Treatment')
    plt.grid(True)
    plt.show()

    # Save the plot
    plt.savefig(f"Result/{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run.png")

