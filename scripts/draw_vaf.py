
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



def get_vaf(repeat_type, cellline):
	df = pd.read_csv("../results/%s_%s_tiers_vaf.txt"%(repeat_type, cellline), sep="\t")
	df_tier2 = df.loc[df["tier2"]=="tier2"]
	df_tier2 = df_tier2.loc[(df_tier2["clipped"] + df_tier2["middle"])>2]
	df_tmp = df_tier2[["vaf"]]
	df_tmp["repeat"] = repeat_type
	df_tmp["cellline"] = cellline
	return df_tmp

repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]


plt.clf()
plt.figure(figsize=(5,3), dpi=1000)
for i in celllines:
	df_tmp = get_vaf("l1", i)
	sns.kdeplot(x=df_tmp["vaf"]*100, label="%s"%i)

plt.xlabel("VAF(%)")
plt.ylabel("Density")
plt.xlim(0,100)
plt.tight_layout()
plt.savefig("vaf_l1.png")


plt.clf()
plt.figure(figsize=(5,3), dpi=1000)
for i in celllines:
	df_tmp = get_vaf("alu", i)
	sns.kdeplot(x=df_tmp["vaf"]*100, label="%s"%i)

plt.xlabel("VAF(%)")
plt.ylabel("Density")
plt.xlim(0,100)
plt.tight_layout()
plt.savefig("vaf_alu.png")


plt.clf()
plt.figure(figsize=(5,3), dpi=1000)
for i in celllines:
	df_tmp = get_vaf("sva", i)
	sns.kdeplot(x=df_tmp["vaf"]*100, label="%s"%i)

plt.xlabel("VAF(%)")
plt.ylabel("Density")
plt.xlim(0,100)
plt.tight_layout()
plt.savefig("vaf_sva.png")


