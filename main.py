from read_my_file import *
def main():
    df_original = dataframe_process("report.pg_matrix 3.tsv")
    df_rerun = dataframe_process("report.pg_matrix-re-run.tsv")
    column_names_rerun = list(df_rerun.columns)[5:-1]
    merged_df = merge_new_run(df_original, df_rerun)
    merged_runs_columns = column_group(merged_df)
    pearson_cor = dataframe_work(merged_df, merged_runs_columns)
    fin_pearson = transform_correlation_dict(pearson_cor)

    # merged_df.to_excel("Manual_result.xlsx", index = False)
    fin_pearson.sort_index().to_excel("Pearson_result.xlsx", index = False)

    experiment = pd.read_excel("With_R1.xlsx")
    experiment_columns = column_group(experiment)

    experiment_pearson = dataframe_work(experiment, experiment_columns)
    experiment_spearman = dataframe_work_spearman(experiment, experiment_columns)

    final_result_pearson = transform_correlation_dict(experiment_pearson)
    final_result_spearman = transform_correlation_dict(experiment_spearman)

    final_result_pearson.sort_values(by='Fraction').to_excel("Final_result_Pearson.xlsx", index = False)
    final_result_spearman.sort_values(by='Fraction').to_excel("Final_result_Spearman.xlsx", index = False)

if __name__ == '__main__':
    main()
