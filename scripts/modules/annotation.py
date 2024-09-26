import os
import subprocess
import pandas as pd


def run_repeatmasker(repeat_type, cellline, caller, df):
    out = open("../tmp/repeatmasker/%s_%s_%s.fa"%(repeat_type, cellline, caller), 'w')
    for i in range(df.shape[0]):
        filename = df.iloc[i]["chrom"] + ";" + df.iloc[i]["pos"] + ";" + df.iloc[i]["caller"]
        out.write(">" + filename + "\n")
        out.write(df.iloc[i]["inserted_sequence"] + "\n")

    out.close()
    command = "RepeatMasker -species human ../tmp/repeatmasker/%s_%s_%s.fa"%(repeat_type, cellline, caller)
    os.system(command)



def run_blast(df):
	blast_result = []
	for i in range(df.shape[0]):
		chrom = df.iloc[i]["chrom"]
		pos = df.iloc[i]["pos"]
		caller = df.iloc[i]["caller"]
		seq = df.iloc[i]["inserted_sequence"]
		name = chrom+"_"+pos+"_"+caller
		out = open("../tmp/blast/tmp.fa", 'w')
		out.write(">%s"%name + "\n")
		out.write(seq + "\n")
		out.close()
		command = "blastn -db ../blast_db/hg38_repeat.fa -query ../tmp/blast/tmp.fa -outfmt 6 "
		result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
		if len(result.stdout) == 0:
			subtype = "."
			match_percent = "."
			overlap_seq = "."
			subject_position = "."
		else:
			df_tmp = pd.DataFrame([row.split("\t") for row in result.stdout.split("\n")])
			subtype = df_tmp.loc[0][1]
			match_percent = df_tmp.loc[0][2]
			overlap_seq = round((float(df_tmp.loc[0][7]) - float(df_tmp.loc[0][6])) / len(seq), 3)
			if float(df_tmp.loc[0][9]) - float(df_tmp.loc[0][8]) > 0:
				overlap_consensus = "+"
			else:
				overlap_consensus = "-"
		output = subtype+";"+str(match_percent)+";"+str(overlap_seq)+";"+subject_position
		blast_result.append(output)
	df["blast"] = blast_result
	return df




