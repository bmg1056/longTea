import os
import pandas as pd

repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG005", "HG00438", "HG02257", "HG02486", "HG02622"]

os.system("mkdir -p ../results/final")

repeat_type = "l1"
cellline = "HG002"

for repeat_type in repeat_types:
    for cellline in celllines:
        df = pd.read_csv("../results/%s_%s_tiers_vaf.txt"%(repeat_type, cellline), sep="\t")
        df = df[["chrom", "pos", "name_dip", "name_sni", "name_mini", "inserted_sequence", "blast(subfamily;match_ratio;insert_start;insert_end;te_start;te_end)", "tsd", "polyA", "direction", "vaf", "tier1", "tier2", "tier3"]]
        df = df.loc[df["tier3"]!="."]
        df = df.rename(columns={"direction":"orientation"})
        df.to_csv("../results/final/%s_%s.txt"%(repeat_type, cellline), sep="\t", index=False)

for repeat_type in repeat_types:
    os.system("cp ../results/%s_hapmapmixture_tier1.txt ../results/final/%s_hapmapmixture_tier1.txt"%(repeat_type, repeat_type))
    os.system("cp ../results/%s_hapmapmixture_tier2.txt ../results/final/%s_hapmapmixture_tier2.txt"%(repeat_type, repeat_type))
    os.system("cp ../results/%s_hapmapmixture_tier3.txt ../results/final/%s_hapmapmixture_tier3.txt"%(repeat_type, repeat_type))




