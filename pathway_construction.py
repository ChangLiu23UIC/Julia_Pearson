from read_my_file import *
def union_lists_to_set(*lists):
    union_set = set()
    for lst in lists:
        union_set.update(lst)
    return union_set

pt_db = pd.read_csv("Top Genes/Top 3 percent Z CCP Genes.csv")

top_ccp_gene_list = list(set(pt_db["Genes"]))





"""
Old runs
"""
# AURKA = df_unique[0]
# AURKB = df_unique[1]
# KIF11 = df_unique[2]
# PRC1 = df_unique[3]
# CCNB1 = df_unique[4]
#
# # Select the proteins based on the KS-SCORE
# # KS_df = pd.read_excel("Sorted_KS-SCORE_with_threshold_0.1_DMSO_vs_whel.xlsx")
# AURKA_df = dmso_shared[dmso_shared["Genes"].isin(AURKA)]
# AURKB_df = dmso_shared[dmso_shared["Genes"].isin(AURKB)]
# KIF11_df = dmso_shared[dmso_shared["Genes"].isin(KIF11)]
# PRC1_df = dmso_shared[dmso_shared["Genes"].isin(PRC1)]
# CCNB1_df = dmso_shared[dmso_shared["Genes"].isin(CCNB1)]
#
# AURKA_df.to_excel("AURKA.xlsx", index=False)
# AURKB_df.to_excel("AURKB.xlsx", index=False)
# KIF11_df.to_excel("KIF11.xlsx", index=False)
# PRC1_df.to_excel("PRC1.xlsx", index=False)
# CCNB1_df.to_excel("CCNB1.xlsx", index=False)
#
# AURKA_list = AURKA_df["Genes"].tolist()
# AURKB_list = AURKB_df["Genes"].tolist()
# KIF11_list = KIF11_df["Genes"].tolist()
# PRC1_list = PRC1_df["Genes"].tolist()
# CCNB1_list = CCNB1_df["Genes"].tolist()
#
# union_set = list(union_lists_to_set(AURKA_list, AURKB_list, KIF11_list, PRC1_list, CCNB1_list))