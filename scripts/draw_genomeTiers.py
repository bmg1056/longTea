
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns

df_easy = pd.read_csv("../genomeTiers/SMaHT_easy_regions_GRCh38_v1.0.bed", sep="\t")
df_difficult = pd.read_csv("../genomeTiers/SMaHT_difficult_regions_all_GRCh38_v1.0.bed", sep="\t")
df_extreme = pd.read_csv("../genomeTiers/SMaHT_extreme_regions_GRCh38_v1.0.bed", sep="\t")

bt_easy = pybedtools.BedTool.from_dataframe(df_easy)
bt_difficult = pybedtools.BedTool.from_dataframe(df_difficult)
bt_extreme = pybedtools.BedTool.from_dataframe(df_extreme)

repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]

def get_overlap(repeat_type, cellline):
	df = pd.read_csv("../results/%s_%s_tiers_vaf.txt"%(repeat_type, cellline), sep="\t")
	df_tier2 = df.loc[df["tier2"]=="tier2"]
	df_tier2 = df_tier2.loc[(df_tier2["clipped"] + df_tier2["middle"])>2]
	df_tier2 = df_tier2[["chrom", "pos", "pos"]]
	bt = pybedtools.BedTool.from_dataframe(df_tier2)
	inter_easy = bt.intersect(bt_easy, wa=True, c=True)
	inter_difficult = bt.intersect(bt_difficult, wa=True, c=True)
	inter_extreme = bt.intersect(bt_extreme, wa=True, c=True)
	df_inter_easy = inter_easy.to_dataframe(names=['chrom', 'start', 'end', 'easy'])
	df_inter_difficult = inter_difficult.to_dataframe(names=['chrom', 'start', 'end', 'difficult'])
	df_inter_extreme = inter_extreme.to_dataframe(names=['chrom', 'start', 'end', 'extreme'])
	df_merge = pd.concat([df_inter_easy, df_inter_difficult[["difficult"]], df_inter_extreme[["extreme"]]], axis=1)
	return [sum(df_merge["easy"]), sum(df_merge["difficult"]), sum(df_merge["extreme"])]

dic_l1 = {}
for i in celllines:
	print(i)
	dic_l1[i] = get_overlap("l1", i)


dic_alu = {}
for i in celllines:
	print(i)
	dic_alu[i] = get_overlap("alu", i)


dic_sva = {}
for i in celllines:
	print(i)
	dic_sva[i] = get_overlap("sva", i)



df_l1 = pd.DataFrame.from_dict(dic_l1, orient="index")
df_alu = pd.DataFrame.from_dict(dic_alu, orient="index")
df_sva = pd.DataFrame.from_dict(dic_sva, orient="index")

df_l1.columns = ["easy", "difficult", "extreme"]
df_alu.columns = ["easy", "difficult", "extreme"]
df_sva.columns = ["easy", "difficult", "extreme"]

df_l1 = df_l1.div(df_l1.sum(axis=1), axis=0)
df_alu = df_alu.div(df_alu.sum(axis=1), axis=0)
df_sva = df_sva.div(df_sva.sum(axis=1), axis=0)

df_l1_melt = pd.melt(df_l1.reset_index(), id_vars=["index"])
df_alu_melt = pd.melt(df_alu.reset_index(), id_vars=["index"])
df_sva_melt = pd.melt(df_sva.reset_index(), id_vars=["index"])


plt.clf()
plt.figure(figsize=(6,4), dpi=1000)
sns.barplot(data=df_l1_melt, x="index", y="value", hue="variable", dodge=False)
plt.xlabel("")
plt.ylabel("Proportion")
plt.tight_layout()
plt.savefig("genomeTiers_l1.png")



plt.clf()
plt.figure(figsize=(6,4), dpi=1000)
sns.barplot(data=df_alu_melt, x="index", y="value", hue="variable", dodge=False)
plt.xlabel("")
plt.ylabel("Proportion")
plt.tight_layout()
plt.savefig("genomeTiers_alu.png")



plt.clf()
plt.figure(figsize=(6,4), dpi=1000)
sns.barplot(data=df_sva_melt, x="index", y="value", hue="variable", dodge=False)
plt.xlabel("")
plt.ylabel("Proportion")
plt.tight_layout()
plt.savefig("genomeTiers_sva.png")





