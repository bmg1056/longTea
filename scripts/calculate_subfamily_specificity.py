import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def tiers(df):
	if len(df["polyA"].split(";")) == 2:
		if df["tsd_original"].split(";")[0] == "tsd":
			tsd_size = int(df["tsd_original"].split(";")[1].split("(")[0])
			if tsd_size > 3 and tsd_size < 31:
				return "tier1"
			else:
				return "tier2"
		else:
			return "tier2"
	else:
		return "."

sub_l1 = ["L1HS", "L1PA2"]
sub_alu = ["AluYa5", "AluYb8", "AluY"]
sub_sva = ["SVA_E", "SVA_F"]

repeat_types = ["l1", "alu", "sva"]
callers = ["dip", "mini", "sni"]

for i in callers:
	df = pd.read_csv("sva_HG002_%s.txt"%i, sep="\t")
	df["sub_rm"] = df["repeatmasker"].str.split(";").str[0]
	df["sub_blast"] = df["blast"].str.split(";").str[0]
	df["inTier"] = df.apply(tiers, axis=1)
	df_rm_specific = df.loc[(df["sub_rm"].isin(sub_sva)) & (~df["sub_blast"].isin(sub_sva))]
	df_blast_specific = df.loc[(~df["sub_rm"].isin(sub_sva)) & (df["sub_blast"].isin(sub_sva))]
	df_overlap = df.loc[(df["sub_rm"].isin(sub_sva)) & (df["sub_blast"].isin(sub_sva))]
	df_rm_specific.shape[0]
	df_rm_specific.loc[df_rm_specific["inTier"]=="tier1"].shape[0]
	df_rm_specific.loc[df_rm_specific["inTier"]=="tier2"].shape[0]
	df_blast_specific.shape[0]
	df_blast_specific.loc[df_blast_specific["inTier"]=="tier1"].shape[0]
	df_blast_specific.loc[df_blast_specific["inTier"]=="tier2"].shape[0]
	df_overlap.shape[0]










i="dip"
df = pd.read_csv("alu_HG002_%s.txt"%i, sep="\t")
df["sub_rm"] = df["repeatmasker"].str.split(";").str[0]
df["sub_blast"] = df["blast"].str.split(";").str[0]
df["inTier"] = df.apply(tiers, axis=1)
df_rm_specific = df.loc[(df["sub_rm"].isin(sub_alu)) & (~df["sub_blast"].isin(sub_alu))]
df_blast_specific = df.loc[(~df["sub_rm"].isin(sub_alu)) & (df["sub_blast"].isin(sub_alu))]
df_overlap = df.loc[(df["sub_rm"].isin(sub_alu)) & (df["sub_blast"].isin(sub_alu))]


## draw repeatmasker score
tmp1 = pd.DataFrame(df_rm_specific["repeatmasker"].str.split(";").str[1].astype(float))
tmp1.columns = ["score"]
tmp1["type"] = "repeatmasker_specific"
tmp2 = pd.DataFrame(df_overlap["repeatmasker"].str.split(";").str[1].astype(float))
tmp2.columns = ["score"]
tmp2["type"] = "overlap"
tmp = pd.concat([tmp1, tmp2], axis=0)

plt.clf()
plt.figure(dpi=1200)
sns.kdeplot(
   data=tmp, x="score", hue="type",
   fill=True, common_norm=False, palette="crest",
   alpha=.5, linewidth=0,
)
plt.tight_layout()
plt.savefig("alu_rm.png")


## draw blast score
tmp1 = pd.DataFrame(df_blast_specific["blast"].str.split(";").str[1].astype(float))
tmp1.columns = ["score"]
tmp1["type"] = "blast_specific"
tmp2 = pd.DataFrame(df_overlap["blast"].str.split(";").str[1].astype(float))
tmp2.columns = ["score"]
tmp2["type"] = "overlap"
tmp = pd.concat([tmp1, tmp2], axis=0)

plt.clf()
plt.figure(dpi=1200)
sns.kdeplot(
   data=tmp, x="score", hue="type",
   fill=True, common_norm=False, palette="crest",
   alpha=.5, linewidth=0,
)
plt.tight_layout()
plt.savefig("alu_blast.png")


## draw length
tmp1 = pd.DataFrame(df_rm_specific.loc[df_rm_specific["inTier"]=="tier1"]["inserted_sequence"].str.len())
tmp1.columns = ["length"]
tmp1["type"] = "repeatmasker_specific"
tmp2 = pd.DataFrame(df_overlap.loc[df_overlap["inTier"]=="tier1"]["inserted_sequence"].str.len())
tmp2.columns = ["length"]
tmp2["type"] = "overlap"
tmp3 = pd.DataFrame(df_blast_specific.loc[df_blast_specific["inTier"]=="tier1"]["inserted_sequence"].str.len())
tmp3.columns = ["length"]
tmp3["type"] = "blast_specific"

tmp = pd.concat([tmp1, tmp2, tmp3], axis=0)

plt.clf()
plt.figure(dpi=1200)
sns.kdeplot(
   data=tmp, x="length", hue="type",
   fill=True, common_norm=False, palette="crest",
   alpha=.5, linewidth=0,
)
plt.xlim(-1000, 2000) 
plt.tight_layout()
plt.savefig("alu_length_tier.png")









