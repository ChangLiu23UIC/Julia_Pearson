from read_my_file import *
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import numpy as np
from pathway_construction import *
from Normalization_methods import *
from scipy.stats import ttest_rel, combine_pvalues
import warnings
from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool, ColumnDataSource, Label



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


# def subset_by_run(df):
#     """
#     Subset my dataframe by run. This looks at the run number and subset each of them into runs.
#     :param df:
#     :return:
#     """
#     genes_column = df[['Genes']]
#
#     x_values = set(col.split('-')[1][1:] for col in df.columns if '-' in col and col != 'Genes')
#
#     subsets = {}
#
#     for x in x_values:
#         # Since there are two types of runs. We need r and n individually.
#         columns_with_x = [col for col in df.columns if f'-r{x}-' in col or f'-n{x}-' in col]
#         columns_with_x_sorted = sorted(columns_with_x, key=lambda col: int(col.split('-F')[-1]))
#         subset = pd.concat([genes_column, df[columns_with_x_sorted]], axis=1)
#         subsets[x] = subset
#     return subsets

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


def ks_test_between_runs(dmso, whel):
    """
    Test the KS score between every runs for DMSO and WHel in a combination. Ex: DMSO run1 vs Whel run1, whelrun2, whelrun3.
    :param dmso:
    :param whel:
    :return:
    """
    genes = dmso["Genes"]
    dmso = dmso.set_index("Genes")
    whel = whel.set_index("Genes")

    dmso_runs = sorted(set(col.split('-')[0] + '-' + col.split('-')[1] for col in dmso.columns))
    whel_runs = sorted(set(col.split('-')[0] + '-' + col.split('-')[1] for col in whel.columns))

    results = []

    for dmso_run in dmso_runs:
        for whel_run in whel_runs:
            for gene in genes:
                dmso_data = dmso.loc[gene, [col for col in dmso.columns if dmso_run in col]].values
                whel_data = whel.loc[gene, [col for col in whel.columns if whel_run in col]].values
                ks_result = ks_2samp(dmso_data, whel_data)
                results.append({
                    'Gene': gene,
                    'Comparison': f"{dmso_run} vs {whel_run}",
                    'KS_statistic': ks_result.statistic,
                    'p_value': ks_result.pvalue
                })

    return pd.DataFrame(results)


def ks_test_total(dmso, whel):
    """
    THe total score of the ks-test into a dataframe between two runs for all genes. THis is for the average scores.
    :param dmso:
    :param whel:
    :return:
    """

    genes = dmso["Genes"]
    ks_results = {'Genes': genes}
    results_statistic = []
    results_pvalue = []

    for gene in genes:
        dmso_data = dmso[dmso["Genes"] == gene].drop(columns=["Genes"]).values.flatten()
        whel_data = whel[whel["Genes"] == gene].drop(columns=["Genes"]).values.flatten()
        ks_result = ks_2samp(dmso_data, whel_data)
        results_statistic.append(ks_result.statistic)
        results_pvalue.append(ks_result.pvalue)

    ks_results['KS_Statistic'] = results_statistic
    ks_results['P_Value'] = results_pvalue
    return pd.DataFrame(ks_results)


def plot_ks_result_histogram(df, method, ks_column='P_Value'):
    """
    Plot the histogram of the ks-score for downstream analysis.
    """
    df['P_Value_First'] = df[ks_column]

    df_sorted = df.sort_values(by='P_Value_First', ascending=False)

    unique_values = df['P_Value_First'].unique()
    num_bins = 27

    plt.figure(figsize=(10, 6))
    plt.hist(df['P_Value_First'], bins=num_bins, edgecolor='black', color='orange')
    plt.title(f'Histogram of the KS test p-values with {method} normalization')
    plt.xlabel('P Value Score')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.savefig(f"hist img/Histogram of the KS test p-values with {method} normalization.png")
    plt.close()

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
    plot the histogram of the data
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

    x_labels = [f"F{i}" for i in range(1,10)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)

    # plt.ylim(-0.2, 1)

    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_mehod} normalized Intensity')
    plt.title(f'{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run')
    plt.legend(title='Treatment')
    plt.grid(True)

    # Save the plot
    plt.savefig(f"img/top z/{normalize_mehod} normalized Intensity for {gene_name} for each DMSO and Whel run.png")
    plt.close()


def plot_diff(DMSO_df, whel_df, gene_name, normalize_method):
    """
    :param dmso:
    :param whel:
    :param gene:
    :param method:
    :return:
    """
    # Filter the DataFrames to only include the specified gene
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

    # Calculate the mean of intensity for each fraction
    sum_df = combined_df.groupby(['Treatment', 'Fraction'])['Intensity'].mean().reset_index()


    # Pivot the DataFrame to have DMSO and whel in separate columns
    pivot_df = sum_df.pivot(index='Fraction', columns='Treatment', values='Intensity').reset_index()


    # Calculate the difference between DMSO and whel
    pivot_df['Difference'] = pivot_df.apply(
        lambda row: ((abs(row['DMSO']) - abs(row['whel'])) / (abs(row['DMSO']) + abs(row["whel"]))) * 100 if  (abs(row['DMSO']) + abs(row["whel"])) != 0 else 0, axis=1
    )

    sns.lineplot(data=pivot_df, x='Fraction', y='Difference', marker='o')
    x_labels = [f"F{i}" for i in range(1,10)]
    plt.xticks(ticks=range(1, 10), labels=x_labels)
    plt.xlabel('Fraction')
    plt.ylabel(f'{normalize_method} normalized Intensity Difference %')
    plt.title(f'{normalize_method} normalized Intensity Difference for {gene_name}')
    plt.grid(True)

    plt.ylim(-100,100)

    # Save the plot
    plt.savefig(f"diff img/{normalize_method} normalized Intensity difference for {gene_name} for each DMSO and Whel run.png")
    plt.close()



def plot_fraction_boxplots(DMSO_df, whel_df, normalize_method):
    """
    Plot box plots of the percentage differences in mean intensities for each fraction between DMSO and whel treatments.

    :param DMSO_df: DataFrame containing DMSO data
    :param whel_df: DataFrame containing whel data
    :param normalize_method: The normalization method used
    :return: None
    """
    # Reshape the DataFrames to long format
    DMSO_long = pd.melt(DMSO_df, id_vars=['Genes'], var_name='Fraction', value_name='Intensity')
    whel_long = pd.melt(whel_df, id_vars=['Genes'], var_name='Fraction', value_name='Intensity')

    # Extract run and fraction information
    DMSO_long[['Treatment', 'Run', 'Fraction']] = DMSO_long['Fraction'].str.extract(r'(DMSO)-n(\d+)-F(\d+)')
    whel_long[['Treatment', 'Run', 'Fraction']] = whel_long['Fraction'].str.extract(r'(whel)-n(\d+)-F(\d+)')

    # Combine the two DataFrames
    combined_df = pd.concat([DMSO_long, whel_long])

    # Convert relevant columns to numeric
    combined_df['Run'] = pd.to_numeric(combined_df['Run'])
    combined_df['Fraction'] = pd.to_numeric(combined_df['Fraction'])

    # Calculate the mean intensities for each gene, fraction, and treatment
    mean_df = combined_df.groupby(['Genes', 'Fraction', 'Treatment'])['Intensity'].median().reset_index()

    # Pivot to get DMSO and whel in separate columns
    pivot_df = mean_df.pivot_table(index=['Genes', 'Fraction'], columns='Treatment', values='Intensity').reset_index()

    # Calculate the percentage difference in mean intensities between DMSO and whel
    pivot_df['Percentage Difference'] = pivot_df.apply(
        lambda row: ((abs(row['DMSO']) - abs(row['whel'])) / (abs(row['DMSO']) + abs(row["whel"]))) * 100 if (abs(
            row['DMSO']) + abs(row["whel"])) != 0 else 0, axis=1
    )

    # Plot the box plots for the percentage differences in mean intensities for each fraction
    plt.figure(figsize=(15, 10))
    sns.lineplot(data=pivot_df, x='Fraction', y='Percentage Difference', markers=True,
                 errorbar ='se', err_style='band')

    plt.ylim(-30, 30)

    plt.xlabel('Fraction', fontsize=15)
    plt.ylabel(f'{normalize_method} normalized Mean Intensity Percentage Difference', fontsize=15)
    plt.title(f'{normalize_method} normalized Mean Intensity Percentage Differences for each Fraction', fontsize=20)
    plt.grid(True)

    # Save the plot
    plt.savefig(f"normalized diff total/{normalize_method}_normalized_Intensity_Percentage_Difference_Boxplot_for_each_Fraction.png")
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


def create_volcano_bokeh_paired(dmso_df, wheldone_df, gene_column='Genes', pseudocount=1e-6):
    # Initialize lists to hold log2 fold changes and p-values
    log2_fc = []
    p_values = []
    genes = []
    fractions = sorted(set(col.split('-')[2] for col in dmso_df.columns if col != gene_column))

    # Find the minimum value in both dataframes to determine the shift needed
    min_value = min(dmso_df.iloc[:, 1:].min().min(), wheldone_df.iloc[:, 1:].min().min())
    shift_value = abs(min_value) + pseudocount if min_value < 0 else pseudocount

    # Create a dictionary to store p-values for each gene
    gene_pvalues = {}

    # Iterate over each gene
    for index, row in dmso_df.iterrows():
        gene = row[gene_column]
        gene_pvalues[gene] = []
        for fraction in fractions:
            # Extract values for the current fraction from both dataframes
            dmso_values = dmso_df.filter(regex=f'DMSO-.*-{fraction}').loc[index].values.astype(float)
            wheldone_values = wheldone_df.filter(regex=f'whel-.*-{fraction}').loc[index].values.astype(float)

            # Add shift value to avoid negative and zero values
            dmso_values = dmso_values + shift_value
            wheldone_values = wheldone_values + shift_value

            # Calculate log2 fold change for this fraction
            log2_fc_value = np.log2(np.mean(wheldone_values) / np.mean(dmso_values))
            log2_fc.append(log2_fc_value)

            # Calculate p-value using paired t-test, handling precision loss warning
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=RuntimeWarning)
                t_stat, p_value = ttest_rel(dmso_values, wheldone_values)

            # Handle NaN p-values
            if np.isnan(p_value):
                p_value = 1.0
            p_values.append(p_value)
            genes.append(gene)
            gene_pvalues[gene].append(p_value)

    # Aggregate p-values for each gene using Fisher's method
    aggregated_pvalues = {}
    for gene, pvals in gene_pvalues.items():
        if pvals:
            _, combined_pvalue = combine_pvalues(pvals)
            aggregated_pvalues[gene] = combined_pvalue

    # Create a dataframe to hold the results
    results_df = pd.DataFrame({
        gene_column: genes,
        'log2_fc': log2_fc,
        'p_value': p_values
    })
    results_df['-log10_p_value'] = -np.log10(results_df['p_value'])

    # Create a dataframe to hold the aggregated results
    aggregated_results = pd.DataFrame({
        gene_column: list(aggregated_pvalues.keys()),
        'combined_p_value': list(aggregated_pvalues.values())
    })
    aggregated_results['-log10_combined_p_value'] = -np.log10(aggregated_results['combined_p_value'])

    # Prepare data for Bokeh
    source = ColumnDataSource(results_df)
    aggregated_source = ColumnDataSource(aggregated_results)

    # Create Bokeh plot
    p = figure(title="Volcano Plot", x_axis_label="Log2 Fold Change", y_axis_label="-Log10 P-value",
               tools="pan,wheel_zoom,box_zoom,reset,save")

    p.scatter('log2_fc', '-log10_p_value', size=8, source=source,
             fill_alpha=0.6, color='grey', legend_label='All genes')

    significant_genes = aggregated_results[aggregated_results['combined_p_value'] < 0.05]
    if not significant_genes.empty:
        source_sig = ColumnDataSource(significant_genes)
        p.scatter('log2_fc', '-log10_combined_p_value', size=8, source=source_sig,
                 fill_alpha=0.6, color='red', legend_label='Significant genes (p < 0.05)')

    hover = HoverTool()
    hover.tooltips = [("Gene", f"@{gene_column}"), ("Log2 Fold Change", "@log2_fc"),
                      ("-Log10 P-value", "@-log10_p_value")]
    p.add_tools(hover)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    output_file("volcano.html")
    show(p)


def create_volcano_t_test(dmso_df, wheldone_df, method, gene_column='Genes', pseudocount=1e-6):
    # Initialize lists to hold log2 fold changes and p-values
    log2_fc = []
    p_values = []
    genes = []

    # Find the minimum value in both dataframes to determine the shift needed
    min_value = min(dmso_df.iloc[:, 1:].min().min(), wheldone_df.iloc[:, 1:].min().min())
    shift_value = abs(min_value) + pseudocount if min_value < 0 else pseudocount

    # Iterate over each gene (each row)
    for index, row in dmso_df.iterrows():
        gene = row[gene_column]
        dmso_values = dmso_df.loc[index].drop(gene_column).values.astype(float)
        wheldone_values = wheldone_df.loc[index].drop(gene_column).values.astype(float)

        # Add shift value to avoid negative and zero values
        dmso_values = dmso_values + shift_value
        wheldone_values = wheldone_values + shift_value

        # Calculate log2 fold change for this gene
        log2_fc_value = np.log2(np.mean(wheldone_values) / np.mean(dmso_values))
        log2_fc.append(log2_fc_value)

        # Calculate p-value using paired t-test, handling precision loss warning
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            t_stat, p_value = ttest_rel(dmso_values, wheldone_values)

        # Handle NaN p-values
        if np.isnan(p_value):
            p_value = 1.0
        p_values.append(p_value)
        genes.append(gene)

    # Create a dataframe to hold the results
    results_df = pd.DataFrame({
        gene_column: genes,
        'log2_fc': log2_fc,
        'p_value': p_values
    })
    results_df['-log10_p_value'] = -np.log10(results_df['p_value'])

    # Prepare data for Bokeh
    source = ColumnDataSource(results_df)

    # Create Bokeh plot
    p = figure(title=f"Volcano Plot for {method}", x_axis_label="Log2 Fold Change (Whel / DMSO)",
               y_axis_label="-Log10 P-value",
               tools="pan,wheel_zoom,box_zoom,reset,save")

    p.scatter('log2_fc', '-log10_p_value', size=8, source=source,
              fill_alpha=0.6, color='grey', legend_label='All genes')

    significant_genes = results_df[results_df['p_value'] < 0.05]
    if not significant_genes.empty:
        source_sig = ColumnDataSource(significant_genes)
        p.scatter('log2_fc', '-log10_p_value', size=8, source=source_sig,
                  fill_alpha=0.6, color='red', legend_label='Significant genes (p < 0.05)')

    hover = HoverTool()
    hover.tooltips = [("Gene", f"@{gene_column}"), ("Log2 Fold Change", "@log2_fc"),
                      ("-Log10 P-value", "@-log10_p_value")]
    p.add_tools(hover)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    output_file(f"Overall_volcano_bokeh/{method}.html")
    show(p)


def create_volcano_bokeh_overall(dmso_df, wheldone_df, ks_df, method, gene_column='Genes', pseudocount=1e-6):
    # Initialize lists to hold log2 fold changes and p-values
    log2_fc = []
    p_values = ks_df["P_Value"].values
    genes = ks_df["Genes"]

    # Find the minimum value in both dataframes to determine the shift needed
    min_value = min(dmso_df.iloc[:, 1:].min().min(), wheldone_df.iloc[:, 1:].min().min())
    shift_value = abs(min_value) + pseudocount if min_value < 0 else pseudocount

    # Iterate over each gene (each row)
    for index, row in dmso_df.iterrows():
        # gene = row[gene_column]
        dmso_values = dmso_df.loc[index].drop(gene_column).values.astype(float)
        wheldone_values = wheldone_df.loc[index].drop(gene_column).values.astype(float)

        # Add shift value to avoid negative and zero values
        dmso_values = dmso_values + shift_value
        wheldone_values = wheldone_values + shift_value

        # Calculate log2 fold change for this gene
        log2_fc_value = np.log2(np.mean(wheldone_values) / np.mean(dmso_values))
        log2_fc.append(log2_fc_value)

    # Check lengths of lists
    print(f"Length of log2_fc: {len(log2_fc)}")
    print(f"Length of p_values: {len(p_values)}")
    print(f"Length of genes: {len(genes)}")

    # Create a dataframe to hold the results
    results_df = pd.DataFrame({
        gene_column: genes,
        'log2_fc': log2_fc,
        'p_value': p_values
    })
    results_df['-log10_p_value'] = -np.log10(results_df['p_value'])

    # Prepare data for Bokeh
    source = ColumnDataSource(results_df)

    # Create Bokeh plot
    p = figure(title=f"Volcano Plot for {method}", x_axis_label="Log2 Fold Change (Whel / DMSO)", y_axis_label="-Log10 P-value of KS-Test",
               tools="pan,wheel_zoom,box_zoom,reset,save")

    p.scatter('log2_fc', '-log10_p_value', size=8, source=source,
             fill_alpha=0.6, color='grey', legend_label='All genes')

    significant_count = (results_df['p_value'] < 0.05).sum()
    significant_genes = results_df[results_df['p_value'] < 0.05]
    if not significant_genes.empty:
        source_sig = ColumnDataSource(significant_genes)
        p.scatter('log2_fc', '-log10_p_value', size=8, source=source_sig,
                 fill_alpha=0.6, color='red', legend_label='Significant genes (p < 0.05)')

    hover = HoverTool()
    hover.tooltips = [("Gene", f"@{gene_column}"), ("Log2 Fold Change", "@log2_fc"),
                      ("-Log10 P-value of KS", "@-log10_p_value")]
    p.add_tools(hover)


    # Add annotation for the number of significant genes
    annotation = Label(x=0.95, y=0.95, x_units='screen', y_units='screen',
                       text=f'Significant genes: {significant_count}')

    p.add_layout(annotation)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    output_file(f"Overall_volcano_bokeh/{method}.html")
    show(p)

def create_ks_volcano_bokeh(ks_df, method, gene_column='Genes'):
    # Initialize lists to hold KS scores and p-values
    ks_scores = 1 / ( 1 - ks_df['KS_Statistic'].values)
    p_values = ks_df['P_Value'].values
    genes = ks_df[gene_column].values

    # Create a dataframe to hold the results
    results_df = pd.DataFrame({
        gene_column: genes,
        'ks_score': ks_scores,
        'p_value': p_values
    })
    results_df['-log10_p_value'] = -np.log10(results_df['p_value'])

    # Prepare data for Bokeh
    source = ColumnDataSource(results_df)

    # Create Bokeh plot
    p = figure(title=f"Volcano Plot (KS Scores) for {method}", x_axis_label="KS Score Modified 1/(1-KS)", y_axis_label="-Log10 P-value",
               tools="pan,wheel_zoom,box_zoom,reset,save")

    p.scatter('ks_score', '-log10_p_value', size=8, source=source,
             fill_alpha=0.6, color='grey', legend_label='All genes')

    significant_genes = results_df[results_df['p_value'] < 0.05]
    if not significant_genes.empty:
        source_sig = ColumnDataSource(significant_genes)
        p.scatter('ks_score', '-log10_p_value', size=8, source=source_sig,
                 fill_alpha=0.6, color='red', legend_label='Significant genes (p < 0.05)')

    hover = HoverTool()
    hover.tooltips = [("Gene", f"@{gene_column}"), ("KS Modified", "@ks_score"),
                      ("-Log10 P-value", "@-log10_p_value")]
    p.add_tools(hover)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    output_file(f"KS_bokeh/{method}_ks.html")
    show(p)


def subset_by_pathway(df, gene_list):
    """
    Subset the dataframe with list, two at once
    :param DMSO:
    :param whel:
    :param gene_list:
    :return:
    """
    sub_df = df[df["Genes"].isin(gene_list)]
    return sub_df


if __name__ == '__main__':
    print("Hello")

    """ 
    # The NSFA normalization
    """
    nsaf_together = nsaf_method(spec_df).dropna().drop_duplicates()
    nsaf_normalized_dmso, nsaf_normalized_whel = separate_dataframe(nsaf_together)

    """
    The TIC normalization
    """
    tic_normalized_dmso = tic_normalization(dmso_shared)
    tic_normalized_whel = tic_normalization(whel_shared)

    # """
    # Median normalization
    # """
    # median_normalized_dmso = median_normalization(dmso_shared)
    # median_normalized_whel = median_normalization(whel_shared)

    """
    Quantile normalization
    """
    quantile_normalized_dmso = quantile_normalization(dmso_shared)
    quantile_normalized_whel = quantile_normalization(whel_shared)

    """
    Z-score normalization
    """
    z_normalized_dmso = z_normalization(dmso_shared)
    z_normalized_whel = z_normalization(whel_shared)

    """
    Variance Stabalize normalization
    """
    var_stab_normalized_dmso = var_stab_normalization(dmso_shared)
    var_stab_normalized_whel = var_stab_normalization(whel_shared)

    """
    Z-Score on NSAF
    """
    nsaf_z_normalized_dmso = z_normalization(nsaf_normalized_dmso)
    nsaf_z_normalized_whel = z_normalization(nsaf_normalized_whel)


    # # Plot normalizations
    hkg = ["AURKA","CCNB1","KIF11","KIFC1","PLK1","PRC1","TACC3", "RIPK1"]
    for i in hkg:
        plot_gene_intensity(nsaf_normalized_dmso, nsaf_normalized_whel, i, "NSAF")
        plot_gene_intensity(tic_normalized_dmso, tic_normalized_whel, i, "TIC")
        plot_gene_intensity(quantile_normalized_dmso, quantile_normalized_whel, i, "Quantile")
        plot_gene_intensity(z_normalized_dmso, z_normalized_whel, i, "Z_norm")
        plot_gene_intensity(var_stab_normalized_dmso, var_stab_normalized_whel, i, "Variance Stabalize")

    for i in hkg:
        plot_diff(nsaf_normalized_dmso, nsaf_normalized_whel, i, 'NSAF')
        plot_diff(quantile_normalized_dmso, quantile_normalized_whel, i, 'Quantile')
        plot_diff(z_normalized_dmso, z_normalized_whel, i, 'Z-Score')
        plot_diff(var_stab_normalized_dmso, var_stab_normalized_whel, i, 'Variance Stabalize')
        plot_diff(tic_normalized_dmso, tic_normalized_whel, i, 'TIC')

    plot_fraction_boxplots(nsaf_normalized_dmso, nsaf_normalized_whel, 'NSAF')
    plot_fraction_boxplots(quantile_normalized_dmso, quantile_normalized_whel, 'Quantile')
    plot_fraction_boxplots(z_normalized_dmso, z_normalized_whel,  'Z-Score')
    plot_fraction_boxplots(var_stab_normalized_dmso, var_stab_normalized_whel, 'Variance Stabalize')
    plot_fraction_boxplots(tic_normalized_dmso, tic_normalized_whel, 'TIC')
    plot_fraction_boxplots(nsaf_z_normalized_dmso, nsaf_z_normalized_whel, "NSAF-Z")
    #
    # plot_fraction_boxplots(subset_by_pathway(nsaf_normalized_dmso, ccp_genes), subset_by_pathway(nsaf_normalized_whel, ccp_genes), 'NSAF_CCP')
    # plot_fraction_boxplots(subset_by_pathway(quantile_normalized_dmso, ccp_genes), subset_by_pathway(quantile_normalized_whel, ccp_genes), 'Quantile_CCP')
    # plot_fraction_boxplots(subset_by_pathway(z_normalized_dmso, ccp_genes), subset_by_pathway(z_normalized_whel, ccp_genes),  'Z-Score_CCP')
    # plot_fraction_boxplots(subset_by_pathway(var_stab_normalized_dmso, ccp_genes), subset_by_pathway(var_stab_normalized_whel, ccp_genes), 'Variance Stabalize_CCP')
    # plot_fraction_boxplots(subset_by_pathway(tic_normalized_dmso, ccp_genes), subset_by_pathway(tic_normalized_whel, ccp_genes), 'TIC_CCP')
    # plot_fraction_boxplots(subset_by_pathway(nsaf_z_normalized_dmso, ccp_genes), subset_by_pathway(nsaf_z_normalized_whel, ccp_genes), "NSAF-Z_CCP")
    # #
    # #
    # ks_z_df = ks_test_total(z_normalized_dmso, z_normalized_whel)
    # ks_quantile_df = ks_test_total(quantile_normalized_dmso, quantile_normalized_whel)
    ks_nsaf_df = ks_test_total(nsaf_normalized_dmso, nsaf_normalized_whel)
    # ks_tic_df = ks_test_total(tic_normalized_dmso, tic_normalized_whel)
    # ks_var_df = ks_test_total(var_stab_normalized_dmso, var_stab_normalized_whel)
    # ks_nsaf_z_df = ks_test_total(nsaf_z_normalized_dmso, nsaf_z_normalized_whel)
    #
    # create_volcano_bokeh_overall(z_normalized_dmso, z_normalized_whel , ks_z_df,"Z")
    # create_volcano_bokeh_overall(quantile_normalized_dmso, quantile_normalized_whel, ks_quantile_df,  "Quantile")
    # create_volcano_bokeh_overall(nsaf_normalized_dmso, nsaf_normalized_whel, ks_nsaf_df, "NSAF")
    # create_volcano_bokeh_overall(tic_normalized_dmso, tic_normalized_whel, ks_tic_df,  "TIC")
    # create_volcano_bokeh_overall(var_stab_normalized_dmso, var_stab_normalized_whel, ks_var_df, "VAR")
    # create_volcano_bokeh_overall(nsaf_z_normalized_dmso, nsaf_z_normalized_whel, ks_nsaf_z_df, "NSAF-Z")
    #
    #
    #
    # #
    # plot_ks_result_histogram(ks_z_df, "Z-Score")
    # plot_ks_result_histogram(ks_quantile_df, "Quantile")
    # plot_ks_result_histogram(ks_nsaf_df, "NSAF")
    # plot_ks_result_histogram(ks_tic_df, "TIC")
    # plot_ks_result_histogram(ks_var_df, "Variance Stabalized")
    # plot_ks_result_histogram(ks_nsaf_z_df, "NSAF-Z")
    # #
    # #
    # # # Get the top 5 percent of the p_value
    top_nsaf = top_5_percent(ks_nsaf_df)
    # top_z = top_5_percent(ks_z_df)
    # top_quantile = top_5_percent(ks_quantile_df)
    # top_tic = top_5_percent(ks_tic_df)
    # top_var = top_5_percent(ks_var_df)
    # top_nsaf_z = top_5_percent(ks_nsaf_z_df)
    #
    # ccp_nsaf = pd.merge(ccp, top_nsaf, on = "Genes", how = "inner")
    # ccp_z = pd.merge(ccp, top_z, on = "Genes", how = "inner")
    # ccp_quantile = pd.merge(ccp, top_quantile, on = "Genes", how = "inner")
    # ccp_tic = pd.merge(ccp, top_tic, on = "Genes", how = "inner")
    # ccp_var = pd.merge(ccp, top_var, on = "Genes", how = "inner")
    # ccp_nsaf_z = pd.merge(ccp, top_nsaf_z, on = "Genes", how = "inner")


    # Pearson correlation of the data
    # Assuming nsaf_normalized_dmso and nsaf_normalized_whel are defined DataFrames
    df1 = nsaf_normalized_dmso.copy()
    df2 = nsaf_normalized_whel.copy()

    # Set the "Genes" column as the index for both DataFrames
    df1.set_index('Genes', inplace=True)
    df2.set_index('Genes', inplace=True)

    # Ensure df2 columns match df1
    df2.columns = df1.columns

    # Aggregate by taking the mean of duplicate entries
    df1 = df1.groupby(df1.index).mean()
    df2 = df2.groupby(df2.index).mean()

    # Initialize a dictionary to store the correlations
    correlations = {}

    # Iterate over each row (gene)
    for gene in df1.index:
        # Select the rows for the gene in both DataFrames
        series1 = df1.loc[gene]
        series2 = df2.loc[gene]

        # Ensure that we are working with Series objects
        if isinstance(series1, pd.DataFrame):
            series1 = series1.squeeze()
        if isinstance(series2, pd.DataFrame):
            series2 = series2.squeeze()

        # Print the types after squeezing
        print(f"Type of series1 after squeezing for gene {gene}: {type(series1)}")
        print(f"Type of series2 after squeezing for gene {gene}: {type(series2)}")

        # Drop NaN values to avoid errors in correlation calculation
        series1 = series1.dropna()
        series2 = series2.dropna()

        # Ensure both series have the same indices after dropping NaNs
        common_idx = series1.index.intersection(series2.index)
        series1 = series1[common_idx]
        series2 = series2[common_idx]

        # Compute Pearson correlation if there are valid indices
        if not series1.empty and not series2.empty:
            correlation = series1.corr(series2)
        else:
            correlation = float('nan')  # Assign NaN if there is insufficient data

        correlations[gene] = correlation

    correlation_df = pd.DataFrame.from_dict(correlations, orient='index', columns=['Pearson Correlation'])
    print(correlation_df)


"""
Seperation
"""
    # create_ks_volcano_bokeh(ks_nsaf_df, "NSAF")
    # create_ks_volcano_bokeh(ks_z_df, "Z")
    # create_ks_volcano_bokeh(ks_quantile_df, "Quantile")
    # create_ks_volcano_bokeh(ks_tic_df, "TIC")
    # create_ks_volcano_bokeh(ks_var_df, "VAR")
    # create_ks_volcano_bokeh(ks_nsaf_z_df, "NSAF-Z")

    # Try to subset the cell cycle pathway proteins o the top 5 percent KS score of normalization
    # top_nsaf_df = top_nsaf[top_nsaf['Genes'].isin(union_set)]
    # top_z_df = top_z[top_z['Genes'].isin(union_set)]
    # top_quantile_df = top_quantile[top_quantile['Genes'].isin(union_set)]
    # top_tic_df = top_tic[top_tic['Genes'].isin(union_set)]
    # top_var_df = top_var[top_var['Genes'].isin(union_set)]


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

