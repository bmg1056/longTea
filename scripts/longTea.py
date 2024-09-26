
import os
import sys
import pandas as pd
import subprocess

import modules.uniform as uniform
import modules.annotation as annotation



repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG005", "HG00438", "HG02257", "HG02486", "HG02622"]
callers = ["dip", "mini", "sni"]

repeat_type = "l1"
cellline = "HG002"
caller = "dip"


os.system("mkdir -p ../tmp")

######## make uniform format
print("-----------------------")
print("1) Make uniform format")
print("\n")


df_dip = uniform.slim_dip(repeat_type, cellline, "dip")
df_mini = uniform.slim_mini(repeat_type, cellline, "mini")
df_sni = uniform.slim_sni(repeat_type, cellline, "sni")


######## run repeatmasker
print("-----------------------")
print("2) Run repeatmasker")
print("\n")


os.system("mkdir -p ../tmp/repeatmasker")

annotation.run_repeatmasker(repeat_type, cellline, "dip", df_dip)
annotation.run_repeatmasker(repeat_type, cellline, "mini", df_mini)
annotation.run_repeatmasker(repeat_type, cellline, "sni", df_sni)


"""
with open("../tmp/repeatmasker/l1_HG002_dip.fa.out") as f:
	tags = f.readlines()

tmp = []
for i in range(3, len(tags)):
	tmp.append([value for value in tags[i].strip().split(" ") if value])

tmp = pd.DataFrame(tmp)
"""


######## run blastn
print("-----------------------")
print("3) Run blastn")
print("\n")


os.system("mkdir -p ../tmp/blast")

df_dip = annotation.run_blast(df_dip)
df_mini = annotation.run_blast(df_mini)
df_sni = annotation.run_blast(df_sni)





