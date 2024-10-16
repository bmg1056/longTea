import pandas as pd


chroms = []
for i in range(1,23):
        chroms.append("chr%s"%i)

chroms.append("chrX"); chroms.append("chrY")


def replace_with_longest(seq):
    sequences = seq.split(",")
    if len(sequences) == 2:
        return max(sequences, key=len)
    else:
        return sequences[0]



def slim_dip(repeat_type, cellline, caller):
    f = open("../calls/%s/%s_%s_%s.vcf"%(repeat_type, cellline, caller, repeat_type), 'r')
    tmp = []
    for k in f:
        if k[0] == "#":
            continue
        else:
            tags = k.strip().split("\t")
            if tags[0] in chroms:
                tmp.append([tags[0], tags[1], "dip", tags[4]])
    f.close()
    df = pd.DataFrame(tmp, columns=["chrom", "pos", "caller", "inserted_sequence"])
    df["inserted_sequence"] = df["inserted_sequence"].apply(replace_with_longest)
    df["name_dip"] = df["chrom"]+"_"+df["pos"].astype(str)+"_"+df["caller"]
    df = df.loc[df.groupby('name_dip')['inserted_sequence'].idxmax()]
    return df


def slim_mini(repeat_type, cellline, caller):
    f = open("../calls/%s/%s_%s_%s.bed"%(repeat_type, cellline, caller, repeat_type), 'r')
    tmp = []
    for k in f:
        if k[0] == "#":
            continue
        else:
            tags = k.strip().split("\t")
            if tags[1] == tags[2]:
                if tags[0] in chroms:
                    tmp.append([tags[0], tags[1], "mini", tags[13]])
    f.close()
    df = pd.DataFrame(tmp, columns=["chrom", "pos", "caller", "inserted_sequence"])
    df["inserted_sequence"] = df["inserted_sequence"].apply(replace_with_longest)
    df["name_mini"] = df["chrom"]+"_"+df["pos"].astype(str)+"_"+df["caller"]
    df = df.loc[df.groupby('name_mini')['inserted_sequence'].idxmax()]
    return df


def slim_sni(repeat_type, cellline, caller):
    f = open("../calls/%s/%s_%s_%s.vcf"%(repeat_type, cellline, caller, repeat_type), 'r')
    tmp = []
    for k in f:
        if k[0] == "#":
            continue
        else:
            tags = k.strip().split("\t")
            if tags[0] in chroms:
                tmp.append([tags[0], tags[1], "sni", tags[4]])
    f.close()
    df = pd.DataFrame(tmp, columns=["chrom", "pos", "caller", "inserted_sequence"])
    df["inserted_sequence"] = df["inserted_sequence"].apply(replace_with_longest)
    df["name_sni"] = df["chrom"]+"_"+df["pos"].astype(str)+"_"+df["caller"]
    df = df.loc[df.groupby('name_sni')['inserted_sequence'].idxmax()]
    return df



def slim_xtea(repeat_type, cellline, caller):
	f = open("../calls/%s/%s_%s_%s.vcf"%(repeat_type, cellline, caller, repeat_type), 'r')
	tmp = []
	for k in f:
		if k[0] == "#":
			continue
		else:
			tags = k.strip().split("\t")
			if tags[0] in chroms:
				tmp.append([tags[0], tags[1], "xtea", tags[4]])
	f.close()
	df = pd.DataFrame(tmp, columns=['chrom', "pos", "caller", "inserted_sequence"])
	df["inserted_sequence"] = df["inserted_sequence"].apply(replace_with_longest)
	df["name_xtea"] = df["chrom"] + "_" + df["pos"].astype(str) + "_" + df["caller"]
	df = df.loc[df.groupby("name_xtea")["inserted_sequence"].idxmax()]
	return df




