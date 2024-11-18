
import pandas as pd




def get_specific(repeat_type):
        df_l1 = pd.read_csv("../results/%s_HG002_tiers_vaf.txt"%repeat_type, sep="\t")
        df_l1_tier1 = df_l1.loc[df_l1["tier1"]=="tier1"]
        df_l1_tier2 = df_l1.loc[df_l1["tier2"]=="tier2"]
        df_l1_tier3 = df_l1.loc[df_l1["tier3"]=="tier3"]
        df_l1_tier1 = df_l1_tier1.loc[(df_l1_tier1["clipped"] + df_l1_tier1["middle"]) > 2]
        df_l1_tier2 = df_l1_tier2.loc[(df_l1_tier2["clipped"] + df_l1_tier2["middle"]) > 2]
        df_l1_tier3 = df_l1_tier3.loc[(df_l1_tier3["clipped"] + df_l1_tier3["middle"]) > 2]
        df_l1_dip = pd.read_csv("../results/%s_HG002_dip.txt"%repeat_type, sep="\t")
        df_l1_mini = pd.read_csv("../results/%s_HG002_mini.txt"%repeat_type, sep="\t")
        df_l1_sni = pd.read_csv("../results/%s_HG002_sni.txt"%repeat_type, sep="\t")
        df_l1_xtea = pd.read_csv("../results/%s_HG002_xtea.txt"%repeat_type, sep="\t")
        df_l1_xtea["subfamily"] = df_l1_xtea["blast"].str.split(";").str[0]
        if repeat_type == "l1":
                df_l1_xtea_tier1 = df_l1_xtea.loc[(df_l1_xtea["subfamily"]=="L1HS") | (df_l1_xtea["subfamily"]=="L1PA2")]
        elif repeat_type == "alu":
                df_l1_xtea_tier1 = df_l1_xtea.loc[(df_l1_xtea["subfamily"]=="AluYa5") | (df_l1_xtea["subfamily"]=="AluYb8") | (df_l1_xtea["subfamily"]=="AluY")]
        if repeat_type == "sva":
                df_l1_xtea_tier1 = df_l1_xtea.loc[(df_l1_xtea["subfamily"]=="SVA_F") | (df_l1_xtea["subfamily"]=="SVA_E")]
        df_tmp = df_l1_tier2[["chrom", "pos"]]
        df_tmp["name"] = df_l1_tier2["chrom"]+";"+df_l1_tier2["pos"].astype(str)
        df_tmp_xtea = df_l1_xtea_tier1[["chrom", "pos", "name_xtea"]]
        df_tmp = df_tmp.sort_values(by="pos")
        df_tmp_xtea = df_tmp_xtea.sort_values(by="pos")
        df_concat = pd.merge_asof(df_tmp, df_tmp_xtea, on="pos", by="chrom", tolerance=30, direction="nearest")
        df_tmp_xtea["is_matched"] = df_tmp_xtea["name_xtea"].isin(df_concat["name_xtea"])
        df_tmp_xtea_unmatched = df_tmp_xtea[df_tmp_xtea["is_matched"]==False].drop(columns=["is_matched"])
        df_concat = pd.concat([df_concat, df_tmp_xtea_unmatched], ignore_index=True)
        df_concat = df_concat.fillna(".")
        df_concat = df_concat.sort_values(by=["chrom", "pos"]).reset_index(drop=True)
        unique_xtea = df_concat.loc[df_concat["name"]=="."]
        unique_cur = df_concat.loc[df_concat["name_xtea"]=="."]
        return unique_cur, unique_xtea

df_l1_longTea, df_l1_xtea = get_specific("l1")
df_alu_longTea, df_alu_xtea = get_specific("alu")
df_sva_longTea, df_sva_xtea = get_specific("sva")



def make_script(repeat, df):
	out = open("script_%s_specific.bat"%repeat, 'w')
	out.write("new" + "\n")
	out.write("snapshotDirectory /home/mib0032/scratch/analysis/smaht_truthsets/figures/%s_specific"%repeat + "\n")
	out.write("genome hg38" + "\n")
	out.write("load /n/storage_test/xiz978/Wangziying/HG002/hifi/HG002_aligned_GRCh38_winnowmap.sorted.bam" + "\n")
	out.write("load /n/storage_test/xiz978/Wangziying/HG002/short_reads/HG002.GRCh38.300x.sorted.bam" + "\n")
	out.write("maxPanelHeight 400" + "\n")
	out.write("\n")
	for i in range(df.shape[0]):
		chrom = df.iloc[i]["chrom"]
		pos = df.iloc[i]["pos"]
		start = int(pos) - 300
		end = int(pos) + 300
		string = chrom+":"+str(start)+"-"+str(end)
		out.write("goto %s"%string + "\n")
		out.write("snapshot %s.png"%string + "\n")
		out.write("\n")
	out.write("exit" + "\n")
	out.close()

make_script("l1_longTea", df_l1_longTea)
make_script("l1_xtea", df_l1_xtea)
make_script("alu_longTea", df_alu_longTea)
make_script("alu_xtea", df_alu_xtea)
make_script("sva_longTea", df_sva_longTea)
make_script("sva_xtea", df_sva_xtea)


