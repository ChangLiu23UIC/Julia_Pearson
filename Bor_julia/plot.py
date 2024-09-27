from main import *

if __name__ == '__main__':
    df_intensity = fill_na_with_half_min(dataframe_process("report.pg_matrix.tsv"))

    df_dmso_nsaf, df_peng_nsaf, df_moxi_nsaf, df_cipr_nsaf = dataframe_process_nsaf("combined_protein.tsv")

    df_dmso_inten, df_peng_inten, df_moxi_inten, df_cipr_inten = seperate_intensity(df_intensity)
    dmso_z, dmso_tic, dmso_med, dmso_quan, dmso_var = all_norm_methods(df_dmso_inten)
    peng_z, peng_tic, peng_med, peng_quan, peng_var = all_norm_methods(df_peng_inten)
    moxi_z, moxi_tic, moxi_med, moxi_quan, moxi_var = all_norm_methods(df_moxi_inten)
    cipr_z, cipr_tic, cipr_med, cipr_quan, cipr_var = all_norm_methods(df_cipr_inten)

    plot_gene_intensity_DP(df_dmso_nsaf, df_peng_nsaf,  "BB_0732", "NSAF")
    plot_gene_intensity_CM(df_dmso_nsaf, df_cipr_nsaf,  "gyrB", "NSAF")
    plot_gene_intensity_DM(df_dmso_nsaf, df_moxi_nsaf,  "gyrB", "NSAF")
    plot_gene_intensity_CM(df_dmso_nsaf, df_cipr_nsaf,  "gyrA", "NSAF")
    plot_gene_intensity_DM(df_dmso_nsaf, df_moxi_nsaf,  "gyrA", "NSAF")

    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplA", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplB", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplC", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rplD", "NSAF")
    # plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsA", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsB", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsC", "NSAF")
    plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsD", "NSAF")


    plot_gene_intensity_DP(dmso_tic, peng_tic,  "BB_0732", "TIC")
    plot_gene_intensity_CM(dmso_tic, cipr_tic,  "gyrB", "TIC")
    plot_gene_intensity_DM(dmso_tic, moxi_tic,  "gyrB", "TIC")
    plot_gene_intensity_CM(dmso_tic, cipr_tic,  "gyrA", "TIC")
    plot_gene_intensity_DM(dmso_tic, moxi_tic,  "gyrA", "TIC")

    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rplA", "TIC")
    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rplB", "TIC")
    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rplC", "TIC")
    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rplD", "TIC")
    # plot_gene_intensity_all(df_dmso_nsaf, df_cipr_nsaf, df_moxi_nsaf, df_peng_nsaf, "rpsA", "NSAF")
    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rpsB", "TIC")
    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rpsC", "TIC")
    plot_gene_intensity_all(dmso_tic, cipr_tic, moxi_tic, peng_tic, "rpsD", "TIC")
