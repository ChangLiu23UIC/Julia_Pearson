from read_my_file_b import *
from normalization_b import *

if __name__ == '__main__':
    df_intensity = fill_na_with_half_min(dataframe_process("report.pg_matrix.tsv"))

    df_dmso_nsaf, df_peng_nsaf, df_moxi_nsaf, df_cipr_nsaf = dataframe_process_nsaf("combined_protein.tsv")

    df_dmso_inten, df_peng_inten, df_moxi_inten, df_cipr_inten = seperate_intensity(df_intensity)
    dmso_z, dmso_tic, dmso_med, dmso_quan, dmso_var = all_norm_methods(df_dmso_inten)
    peng_z, peng_tic, peng_med, peng_quan, peng_var = all_norm_methods(df_peng_inten)
    moxi_z, moxi_tic, moxi_med, moxi_quan, moxi_var = all_norm_methods(df_moxi_inten)
    cipr_z, cipr_tic, cipr_med, cipr_quan, cipr_var = all_norm_methods(df_dmso_inten)