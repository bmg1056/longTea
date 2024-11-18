
import pandas as pd

celllines = ["HG005", "HG02622", "HG02486", "HG02257", "HG002", "HG00438"]
repeat_type = ["l1", "alu", "sva"]
tiers = ["tier1", "tier2", "tier3"]


def merge_two_files(df1_tmp, df2, name2):
    df2_tmp = df2[["chrom", "pos"]]
    df2_tmp["%s"%name2] = df2_tmp["chrom"]+":"+df2_tmp["pos"].astype(str)
    df1_tmp = df1_tmp.sort_values(by="pos")
    df2_tmp = df2_tmp.sort_values(by="pos")
    df_concat = pd.merge_asof(df1_tmp, df2_tmp, on="pos", by="chrom", tolerance=30, direction="nearest")
    df2_tmp["is_matched"] = df2_tmp["%s"%name2].isin(df_concat["%s"%name2])
    df2_tmp_unmatched = df2_tmp[df2_tmp["is_matched"]== False].drop(columns=["is_matched"])
    df_concat = pd.concat([df_concat, df2_tmp_unmatched], ignore_index=True)
    return df_concat

def make_hapmapmixture_truthset(df_dic, tier):
    df_concat = df_dic["HG005"][tier][["chrom", "pos"]]
    df_concat["HG005"] = df_concat["chrom"]+":"+df_concat["pos"].astype(str)
    df_concat = merge_two_files(df_concat, df_dic["HG02622"][tier], "HG02622")
    df_concat = merge_two_files(df_concat, df_dic["HG02486"][tier], "HG02486")
    df_concat = merge_two_files(df_concat, df_dic["HG02257"][tier], "HG02257")
    df_concat = merge_two_files(df_concat, df_dic["HG002"][tier], "HG002")
    df_concat = merge_two_files(df_concat, df_dic["HG00438"][tier], "HG00438")
    df_concat = df_concat.fillna(".")
    df_concat = df_concat.sort_values(by=["chrom", "pos"]).reset_index(drop=True)
    return df_concat

def insertion_type(x):
    if x["HG005"] == ".":
        return "somatic"
    else:
        return "germline"

def make_dictionary(repeat_type):
    df_dic = {}
    for i in celllines:
        tmp = pd.read_csv("../results/%s_%s_tiers.txt"%(repeat_type, i), sep="\t")
        tmp["name_cellline"] = tmp["chrom"]+";"+tmp["pos"].astype("str")+";"+"%s"%i
        df_dic[i] = {}
        for j in tiers:
            tmp1 = tmp.loc[tmp[j]==j]
            df_dic[i][j] = tmp1
    return df_dic


for i in repeat_type:
    df_dic = make_dictionary(i)
    df_tier1 = make_hapmapmixture_truthset(df_dic, "tier1")
    df_tier2 = make_hapmapmixture_truthset(df_dic, "tier2")
    df_tier3 = make_hapmapmixture_truthset(df_dic, "tier3")
    df_tier1["insertion_type"] = df_tier1.apply(insertion_type, axis=1)
    df_tier2["insertion_type"] = df_tier2.apply(insertion_type, axis=1)
    df_tier3["insertion_type"] = df_tier3.apply(insertion_type, axis=1)
    df_tier1.to_csv("../results/%s_hapmapmixture_tier1.txt"%(i), sep="\t", index=False)
    df_tier2.to_csv("../results/%s_hapmapmixture_tier2.txt"%(i), sep="\t", index=False)
    df_tier3.to_csv("../results/%s_hapmapmixture_tier3.txt"%(i), sep="\t", index=False)



