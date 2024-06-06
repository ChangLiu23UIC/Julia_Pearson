import pandas as pd


def unique_columns(df):
    return [list(set(df[col])) for col in df.columns]

def union_lists_to_set(*lists):
    union_set = set()
    for lst in lists:
        union_set.update(lst)
    return union_set

pt_db = pd.read_excel("pathway-proteins_all.xlsx")

df_unique = unique_columns(pt_db)

AURKA = df_unique[0]
AURKB = df_unique[1]
KIF11 = df_unique[2]
PRC1 = df_unique[3]
CCNB1 = df_unique[4]

# union_set = union_lists_to_set(AURKA, AURKB, KIF11, PRC1, CCNB1)
