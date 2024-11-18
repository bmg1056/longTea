
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]

def get_number(repeat_type, cellline):
	df = pd.read_csv("../results/%s_%s_tiers_vaf.txt"%(repeat_type, cellline), sep="\t")
	df_tier1 = df.loc[df["tier1"]=="tier1"]
	df_tier2 = df.loc[df["tier2"]=="tier2"]
	df_tier3 = df.loc[df["tier3"]=="tier3"]
	df_tier1 = df_tier1.loc[(df_tier1["clipped"] + df_tier1["middle"])>2]
	df_tier2 = df_tier2.loc[(df_tier2["clipped"] + df_tier2["middle"])>2]
	df_tier3 = df_tier3.loc[(df_tier3["clipped"] + df_tier3["middle"])>2]
	return [df_tier1.shape[0], df_tier2.shape[0], df_tier3.shape[0]]



dic_l1 = {}
for i in celllines:
	print(i)
	dic_l1[i] = get_number("l1", i)


dic_alu = {}
for i in celllines:
	print(i)
	dic_alu[i] = get_number("alu", i)


dic_sva = {}
for i in celllines:
	print(i)
	dic_sva[i] = get_number("sva", i)



df_l1 = pd.DataFrame.from_dict(dic_l1, orient="index")
df_alu = pd.DataFrame.from_dict(dic_alu, orient="index")
df_sva = pd.DataFrame.from_dict(dic_sva, orient="index")

df_l1.columns = ["tier1", "tier2", "tier3"]
df_alu.columns = ["tier1", "tier2", "tier3"]
df_sva.columns = ["tier1", "tier2", "tier3"]

df_l1_melt = pd.melt(df_l1.reset_index(), id_vars=["index"])
df_alu_melt = pd.melt(df_alu.reset_index(), id_vars=["index"])
df_sva_melt = pd.melt(df_sva.reset_index(), id_vars=["index"])


plt.clf()
plt.figure(figsize=(6,5), dpi=1000)
sns.barplot(data=df_l1_melt, x="index", y="value", hue="variable")
plt.xlabel("")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("count_l1.png")



plt.clf()
plt.figure(figsize=(6,5), dpi=1000)
sns.barplot(data=df_alu_melt, x="index", y="value", hue="variable")
plt.xlabel("")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("count_alu.png")



plt.clf()
plt.figure(figsize=(6,5), dpi=1000)
sns.barplot(data=df_sva_melt, x="index", y="value", hue="variable")
plt.xlabel("")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("count_sva.png")



