
import os
import sys
import pandas as pd
import subprocess

import modules.uniform as uniform
import modules.repeatmasker as repeatmasker



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

repeatmasker.run_repeatmasker(repeat_type, cellline, "dip", df_dip)
repeatmasker.run_repeatmasker(repeat_type, cellline, "mini", df_mini)
repeatmasker.run_repeatmasker(repeat_type, cellline, "sni", df_sni)







######## run blastn
print("-----------------------")
print("3) Run blastn")
print("\n")







out = open("tmp.fa", 'w')
out.write(">%s"%"test" + "\n")
out.write(df_dip.iloc[0]["inserted_sequence"] + "\n")
out.close()

command = "blastn -db ../blast_db/hg38_repeat.fa -query tmp.fa -outfmt 6 "
result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
result.stdout


os.system("blastn -db ./blast_db/hg38_repeat.fa -query ./blast_tmp/%s/%s_%s/%s.fa -outfmt 6 -out ./blast_tmp/%s/%s_%s/%s.blast"%(repeat_type, cellline, caller, filename, repeat_type, cellline, caller, filename))


