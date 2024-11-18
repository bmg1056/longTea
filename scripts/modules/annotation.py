import os
import subprocess
import pandas as pd

from Bio import SeqIO

import re


def run_repeatmasker(repeat_type, cellline, caller, df):
    out = open("../tmp/repeatmasker/%s_%s_%s.fa"%(repeat_type, cellline, caller), 'w')
    for i in range(df.shape[0]):
        filename = df.iloc[i]["chrom"] + "_" + df.iloc[i]["pos"] + "_" + df.iloc[i]["caller"]
        out.write(">" + filename + "\n")
        out.write(df.iloc[i]["inserted_sequence"] + "\n")

    out.close()
    command = "RepeatMasker -species human ../tmp/repeatmasker/%s_%s_%s.fa"%(repeat_type, cellline, caller)
    os.system(command)



def add_repeatmasker(repeat_type, cellline, caller, df):
        with open("../tmp/repeatmasker/%s_%s_%s.fa.out"%(repeat_type, cellline, caller)) as f:
                tags = f.readlines()
        tmp = []
        for i in range(3, len(tags)):
                tmp.append([value for value in tags[i].strip().split(" ") if value])
        tmp = pd.DataFrame(tmp)
        dic = {}
        for i in set(tmp[4]):
                if repeat_type == "l1":
                        df_tmp = tmp.loc[(tmp[4]==i) & (tmp[10]=="LINE/L1")]
                elif repeat_type == "alu":
                        df_tmp = tmp.loc[(tmp[4]==i) & (tmp[10]=="SINE/Alu")]
                elif repeat_type == "sva":
                        df_tmp = tmp.loc[(tmp[4]==i) & (tmp[10]=="Retroposon/SVA")]
                if len(df_tmp[9]) == 0:
                        dic[i] = ".;.;.;.;.;."
                else:
                		if df_tmp[11].values[0][0] == "(":
		                        dic[i] = df_tmp[9].values[0]+";"+str(df_tmp[0].values[0])+";"+str(df_tmp[5].values[0])+";"+str(df_tmp[6].values[0])+";"+str(df_tmp[12].values[0])+";"+str(df_tmp[13].values[0])
                		else:
		                        dic[i] = df_tmp[9].values[0]+";"+str(df_tmp[0].values[0])+";"+str(df_tmp[5].values[0])+";"+str(df_tmp[6].values[0])+";"+str(df_tmp[11].values[0])+";"+str(df_tmp[12].values[0])
        results_repeatmasker = []
        for i in range(df.shape[0]):
                if df.iloc[i]["name_%s"%caller] in dic.keys():
                        results_repeatmasker.append(dic[df.iloc[i]["name_%s"%caller]])
                else:
                        results_repeatmasker.append(".")
        df["repeatmasker"] = results_repeatmasker
        return df







def run_blast(df, repeat_type, cellline):
        blast_result = []
        for i in range(df.shape[0]):
                chrom = df.iloc[i]["chrom"]
                pos = df.iloc[i]["pos"]
                caller = df.iloc[i]["caller"]
                seq = df.iloc[i]["inserted_sequence"]
                name = chrom+"_"+pos+"_"+caller
                out = open("../tmp/blast/tmp_%s_%s.fa"%(repeat_type, cellline), 'w')
                out.write(">%s"%name + "\n")
                out.write(seq + "\n")
                out.close()
                command = "blastn -db ../blast_db/hg38_repeat.fa -query ../tmp/blast/tmp_%s_%s.fa -outfmt 6 "%(repeat_type, cellline)
                result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
                if len(result.stdout) == 0:
                        subtype = "."
                        bit_score = "."
                        overlap_insert_start = "."
                        overlap_insert_end = "."
                        overlap_te_start = "."
                        overlap_te_end = "."
                else:
                        df_tmp = pd.DataFrame([row.split("\t") for row in result.stdout.split("\n")])
                        subtype = df_tmp.loc[0][1]
                        bit_score = df_tmp.loc[0][11]
                        overlap_insert_start = str(df_tmp.loc[0][6])
                        overlap_insert_end = str(df_tmp.loc[0][7])
                        overlap_te_start = str(df_tmp.loc[0][8])
                        overlap_te_end = str(df_tmp.loc[0][9])
                output = subtype+";"+str(bit_score)+";"+str(overlap_insert_start)+";"+str(overlap_insert_end)+";"+str(overlap_te_start)+";"+str(overlap_te_end)
                blast_result.append(output)
        df["blast"] = blast_result
        return df







def get_cigar(df, repeat_type, cellline):
        ref = list(SeqIO.parse("/home/ch252274/work/refs/gatk_bundle/hg38/Homo_sapiens_assembly38.fasta", "fasta"))

        list_cigar = []
        for i in range(df.shape[0]):
                chrom = df.iloc[i]["chrom"][3:]
                pos = int(df.iloc[i]["pos"])
                seq = df.iloc[i]["inserted_sequence"][1:]
                if chrom == "X":
                        chrom = 22
                elif chrom == "Y":
                        chrom = 23
                else:
                        chrom = int(chrom) - 1

                ref_left = ref[chrom].seq[pos-500:pos]
                ref_right = ref[chrom].seq[pos:pos+500]

                insertion_left = seq[:100]
                insertion_right = seq[-100:]
                pseudoread_left = ref_left + insertion_left
                pseudoread_right = insertion_right + ref_right
                ref_seq = ref_left + ref_right

                out_ref = open("../tmp/bwa/tmp_%s_%s.fa"%(repeat_type, cellline), 'w')
                out_ref.write(">tmp" + "\n")
                out_ref.write(str(ref_seq))
                out_ref.close()
                os.system("bwa index ../tmp/bwa/tmp_%s_%s.fa"%(repeat_type, cellline))

                out_read = open("../tmp/bwa/tmp_%s_%s.fq"%(repeat_type, cellline), 'w')
                out_read.write(">left" + "\n")
                out_read.write(str(pseudoread_left) + "\n")
                out_read.write("+" + "\n")
                out_read.write("I" * len(pseudoread_left) + "\n")
                out_read.write(">right" + "\n")
                out_read.write(str(pseudoread_right) + "\n")
                out_read.write("+" + "\n")
                out_read.write("I" * len(pseudoread_right) + "\n")
                out_read.close()

                command = "bwa mem -v 1 ../tmp/bwa/tmp_%s_%s.fa ../tmp/bwa/tmp_%s_%s.fq"%(repeat_type, cellline, repeat_type, cellline)
                result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
                df_tmp = pd.DataFrame([row.split("\t") for row in result.stdout.split("\n")])
                cigar_left = ";".join(df_tmp.loc[df_tmp[0] == "left"][5].values)
                cigar_right = ";".join(df_tmp.loc[df_tmp[0] == "right"][5].values)
                list_cigar.append(cigar_left+"/"+cigar_right)
        df["cigar"] = list_cigar
        return df



################ TSD


def process_string(s):
#    matches = re.findall(r'(\d+)M(?:\d{1}+[DI])?(\d+)?M', s)
    matches = re.findall(r'(\d+)M(?:\d{1}[DI])?(\d+)?M', s)
    if matches:
        # Sum the values of 'M' parts
        m_total = sum(int(m) for m in matches[0] if m)
        # Find the S value
        s_match = re.search(r'(\d+)S', s)
        s_value = s_match.group(1) if s_match else '0'
        return f'{m_total}M{s_value}S'
    return s




def get_tsd(df):
    pattern_left = r'^(\d+M\d+S)$|(\d+M\d+[ID]\d+M\d+S)'
    pattern_right = r'^(\d+S\d+M)$|(\d+S\d+M\d+[ID]\d+M)'
    cigars_left = [s for s in df["cigar"].split("/")[0].split(";") if re.match(pattern_left, s)]
    cigars_right = [s for s in df["cigar"].split("/")[1].split(";") if re.match(pattern_right, s)]
    if len(cigars_left) == 1 and len(cigars_right) == 1:
        cigar_left = re.findall(r'\d+[A-Z]', cigars_left[0])
        cigar_right = re.findall(r'\d+[A-Z]', cigars_right[0])
        cigar_left_size = sum(int(item[:-1]) for item in cigar_left if item.endswith('M'))
        cigar_right_size = sum(int(item[:-1]) for item in cigar_right if item.endswith('M'))
#        tsd_size = abs(float(cigar_right[-1][:-1]) - float(cigar_left[0][:-1]))
        if cigar_right_size > cigar_left_size:
            tsd_size = float(cigar_right_size) - 500
        else:
            tsd_size = float(cigar_left_size) - 500
        if float(cigar_right_size) <= float(cigar_left_size):
            tsd_seq = df["inserted_sequence"][:int(tsd_size)]
        else:
            tsd_seq = df["inserted_sequence"][len(df["inserted_sequence"])-int(tsd_size):len(df["inserted_sequence"])]
        return "tsd;%s(%s)"%(str(int(tsd_size)), tsd_seq)
    else:
        return "no_tsd;no_appropriate_cigar(%s|%s)"%(str(len(cigars_left)), str(len(cigars_right)))





def get_tsd_original(df):
    pattern_left = r'^(\d+M\d+S)'
    pattern_right = r'^(\d+S\d+M)'
    cigars_left = [s for s in df["cigar"].split("/")[0].split(";") if re.match(pattern_left, s)]
    cigars_right = [s for s in df["cigar"].split("/")[1].split(";") if re.match(pattern_right, s)]
    if len(cigars_left) == 1 and len(cigars_right) == 1:
        cigar_left = re.findall(r'\d+[A-Z]', cigars_left[0])
        cigar_right = re.findall(r'\d+[A-Z]', cigars_right[0])
#        tsd_size = abs(float(cigar_right[-1][:-1]) - float(cigar_left[0][:-1]))
        if float(cigar_right[-1][:-1]) > float(cigar_left[0][:-1]):
            tsd_size = float(cigar_right[-1][:-1]) - 500
        else:
            tsd_size = float(cigar_left[0][:-1]) - 500
        if float(cigar_right[-1][:-1]) <= float(cigar_left[0][:-1]):
            tsd_seq = df["inserted_sequence"][:int(tsd_size)]
        else:
            tsd_seq = df["inserted_sequence"][len(df["inserted_sequence"])-int(tsd_size):len(df["inserted_sequence"])]
        return "tsd;%s(%s)"%(str(int(tsd_size)), tsd_seq)
    else:
        return "no_tsd;no_appropriate_cigar(%s|%s)"%(str(len(cigars_left)), str(len(cigars_right)))





################## polyA





def detect_polyA_signal(sequence, length=10, allowed_mismatches=2):
    seq_left = sequence.split(";")[0]
    seq_right = sequence.split(";")[1]
    polyA_positions = []
    polyT_positions = []
    for i in range(len(seq_left) - length + 1):
        subseq = seq_left[i:i + length]
        mismatch_count = sum(1 for base in subseq if base != 'T')
        if mismatch_count <= allowed_mismatches:
            polyT_positions.append([i, subseq])
    for j in range(len(seq_right) - length + 1):
        subseq = seq_right[j:j + length]
        mismatch_count = sum(1 for base in subseq if base != 'A')
        if mismatch_count <= allowed_mismatches:
            polyA_positions.append([j, subseq])
    if len(polyA_positions) > 0 and len(polyT_positions) > 0:
        if polyA_positions[0][0] < polyT_positions[0][0]:
            return "polyA;%s"%(polyA_positions[0][1][::-1])
        else:
            return "polyT;%s"%(polyT_positions[0][1])
    elif len(polyA_positions) > 0:
        return "polyA;%s"%(polyA_positions[0][1][::-1])
    elif len(polyT_positions) > 0:
        return "polyT;%s"%(polyT_positions[0][1])
    else:
        return "no_polyAT"


def get_polyA(df):
    if df["tsd_original"].split(";")[0] == "tsd":
        pattern_left = r'^(\d+M\d+S)'
        pattern_right = r'^(\d+S\d+M)'
        cigars_left = [s for s in df["cigar"].split("/")[0].split(";") if re.match(pattern_left, s)]
        cigars_right = [s for s in df["cigar"].split("/")[1].split(";") if re.match(pattern_right, s)]
        cigar_left = re.findall(r'\d+[A-Z]', cigars_left[0])
        cigar_right = re.findall(r'\d+[A-Z]', cigars_right[0])
        if float(cigar_right[-1][:-1]) <= float(cigar_left[0][:-1]):
            size = int(cigar_left[0][:-1]) - 500
            seq = df["inserted_sequence"][size:size+20]+";"+df["inserted_sequence"][-20:][::-1]
        else:
            size = int(cigar_right[-1][:-1]) - 500
            seq = df["inserted_sequence"][0:20]+";"+df["inserted_sequence"][-20-size:-size][::-1]
    else:
        seq = df["inserted_sequence"][0:20]+";"+df["inserted_sequence"][-20:][::-1]
    return detect_polyA_signal(seq)







