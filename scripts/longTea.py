
import os
import sys
import pandas as pd
import subprocess

import modules.uniform as uniform
import modules.annotation as annotation
import modules.concat as concat






def run_longTea(repeat_type, cellline):
        os.system("mkdir -p ../tmp")
        ######## make uniform format
        print("-----------------------")
        print("1) Make uniform format (%s, %s)"%(repeat_type, cellline))
        print("\n")

        df_dip = uniform.slim_dip(repeat_type, cellline, "dip")
        df_mini = uniform.slim_mini(repeat_type, cellline, "mini")
        df_sni = uniform.slim_sni(repeat_type, cellline, "sni")

######## Add repeatmasker
#print("-----------------------")
#print("2) Add repeatmasker")
#print("\n")


#os.system("mkdir -p ../tmp/repeatmasker")

#annotation.run_repeatmasker(repeat_type, cellline, "dip", df_dip)
#annotation.run_repeatmasker(repeat_type, cellline, "mini", df_mini)
#annotation.run_repeatmasker(repeat_type, cellline, "sni", df_sni)


#with open("../tmp/repeatmasker/l1_HG002_dip.fa.out") as f:
#        tags = f.readlines()

#tmp = []
#for i in range(3, len(tags)):
#        tmp.append([value for value in tags[i].strip().split(" ") if value])

#tmp = pd.DataFrame(tmp)

        ######## Add blastn
        print("-----------------------")
        print("3) Add blastn (%s, %s)"%(repeat_type, cellline))
        print("\n")
        os.system("mkdir -p ../tmp/blast")
        df_dip = annotation.run_blast(df_dip)
        df_mini = annotation.run_blast(df_mini)
        df_sni = annotation.run_blast(df_sni)

        ######## Add cigar
        print("-----------------------")
        print("4) Add cigar (%s, %s)"%(repeat_type, cellline))
        print("\n")
        os.system("mkdir -p ../tmp/bwa")
        df_dip = annotation.get_cigar(df_dip)
        df_mini = annotation.get_cigar(df_mini)
        df_sni = annotation.get_cigar(df_sni)

        ######## Get TPRT signal
        print("-----------------------")
        print("5) Get TPRT signal (%s, %s)"%(repeat_type, cellline))
        print("\n")
        df_dip["tsd"] = df_dip.apply(annotation.get_tsd, axis=1)
        df_mini["tsd"] = df_mini.apply(annotation.get_tsd, axis=1)
        df_sni["tsd"] = df_sni.apply(annotation.get_tsd, axis=1)
        df_dip["polyA"] = df_dip.apply(annotation.get_polyA, axis=1)
        df_mini["polyA"] = df_mini.apply(annotation.get_polyA, axis=1)
        df_sni["polyA"] = df_sni.apply(annotation.get_polyA, axis=1)

        ######## save dataframe
        print("-----------------------")
        print("6) Save data (%s, %s)"%(repeat_type, cellline))
        print("\n")
        os.system("mkdir -p ../results")
        df_dip.to_csv("../results/%s_%s_dip.txt"%(repeat_type, cellline), sep="\t", index=False)
        df_mini.to_csv("../results/%s_%s_mini.txt"%(repeat_type, cellline), sep="\t", index=False)
        df_sni.to_csv("../results/%s_%s_sni.txt"%(repeat_type, cellline), sep="\t", index=False)

        ######## concat results
        print("-----------------------")
        print("7) Concat results (%s, %s)"%(repeat_type, cellline))
        print("\n")
        df_dip = pd.read_csv("../results/%s_%s_dip.txt"%(repeat_type, cellline), sep="\t")
        df_mini = pd.read_csv("../results/%s_%s_mini.txt"%(repeat_type, cellline), sep="\t")
        df_sni = pd.read_csv("../results/%s_%s_sni.txt"%(repeat_type, cellline), sep="\t")
        df = concat.concat_files(df_dip, df_mini, df_sni, repeat_type, cellline)
        df.to_csv("../results/%s_%s_tiers.txt"%(repeat_type, cellline), sep="\t", index=False)





def run_longTea_xtea(repeat_type, cellline):
        os.system("mkdir -p ../tmp")
        ######## make uniform format
        print("-----------------------")
        print("1) Make uniform format (%s, %s)"%(repeat_type, cellline))
        print("\n")

        df_xtea = uniform.slim_xtea(repeat_type, cellline, "xtea")


        ######## Add blastn
        print("-----------------------")
        print("3) Add blastn (%s, %s)"%(repeat_type, cellline))
        print("\n")
        os.system("mkdir -p ../tmp/blast")
        df_xtea = annotation.run_blast(df_xtea)

        ######## Add cigar
        print("-----------------------")
        print("4) Add cigar (%s, %s)"%(repeat_type, cellline))
        print("\n")
        os.system("mkdir -p ../tmp/bwa")
        df_xtea = annotation.get_cigar(df_xtea)

        ######## Get TPRT signal
        print("-----------------------")
        print("5) Get TPRT signal (%s, %s)"%(repeat_type, cellline))
        print("\n")
        df_xtea["tsd"] = df_xtea.apply(annotation.get_tsd, axis=1)
        df_xtea["polyA"] = df_xtea.apply(annotation.get_polyA, axis=1)

        ######## save dataframe
        print("-----------------------")
        print("6) Save data (%s, %s)"%(repeat_type, cellline))
        print("\n")
        os.system("mkdir -p ../results")
        df_xtea.to_csv("../results/%s_%s_xtea.txt"%(repeat_type, cellline), sep="\t", index=False)




repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG005", "HG00438", "HG02257", "HG02486", "HG02622"]

for i in celllines:
        for j in repeat_types:
            run_longTea(j, i)


for i in repeat_types:
    run_longTea_xtea(i, "HG002")


