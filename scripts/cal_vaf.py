import os
import pandas as pd
import pysam
from collections import Counter


def get_vaf(name, df_caller, cellline):
	samfile = pysam.AlignmentFile("/home/mib0032/scratch/data/pangenome/%s_hifi.bam"%cellline, "rb")
	chrom = name.split("_")[0]
	pos = int(name.split("_")[1])
	start = pos - 100
	end = pos + 100
	inserted_sequence =  df_caller.loc[name]["inserted_sequence"]
	length_inserted_sequence = len(inserted_sequence)
	overlapped_reads = samfile.fetch(chrom, start, end)
	# extract the reads
	dic_reads = {}
	for read in overlapped_reads:
		if not read.is_secondary and read.mapq>=20:
			if read.query_name in dic_reads.keys():
				read_hardclip = 0
				original_hardclip = 0
				if read.cigar[0][0] == 5 or read.cigar[-1][0] == 5:
					read_hardclip = 1
				if dic_reads[read.query_name].cigar[0][0] == 5 or dic_reads[read.query_name].cigar[-1][0] == 5:
					original_hardclip = 1
				if read_hardclip == 1 and original_hardclip == 0:
					continue
				elif read_hardclip == 0 and original_hardclip == 1:
					dic_reads[read.query_name] = read
				elif read_hardclip == 1 and original_hardclip == 1:
					if int(read.mapq) > int(dic_reads[read.query_name].mapq):
						dic_reads[read.query_name] = read
					else:
						continue
				else:
					if int(read.mapq) > int(dic_reads[read.query_name].mapq):
						dic_reads[read.query_name] = read
					else:
						continue			
			else:
				dic_reads[read.query_name] = read
	# extract the reads spanning the breakpoint
	remove_list = []
	for i in dic_reads.keys():
		start_pos = dic_reads[i].reference_start
		end_pos = dic_reads[i].reference_end
		if dic_reads[i].cigartuples[0][0] == 4:
			start_pos -= dic_reads[i].cigartuples[0][1]
		if dic_reads[i].cigartuples[-1][0] == 4:
			end_pos += dic_reads[i].cigartuples[-1][1]
		if start_pos <= pos and pos <= end_pos:
			continue
		else:
			remove_list.append(i)
	dic_reads = {key: value for key, value in dic_reads.items() if key not in remove_list}
	# calculate the boundary
	left_boundary = []
	right_boundary = []
	for i in dic_reads.keys():
		if dic_reads[i].cigar[0][0] == 4 and dic_reads[i].cigar[-1][0] == 4 and dic_reads[i].cigar[0][1] >= 10 and dic_reads[i].cigar[-1][1] < 10:
			left_boundary.append(dic_reads[i].reference_start)
		elif dic_reads[i].cigar[0][0] == 4 and dic_reads[i].cigar[-1][0] != 4 and dic_reads[i].cigar[0][1] >= 10:
			left_boundary.append(dic_reads[i].reference_start)
		elif dic_reads[i].cigar[0][0] == 4 and dic_reads[i].cigar[-1][0] == 4 and dic_reads[i].cigar[0][1] < 10 and dic_reads[i].cigar[-1][1] >= 10:
			right_boundary.append(dic_reads[i].reference_end)
		elif dic_reads[i].cigar[-1][0] == 4 and dic_reads[i].cigar[0][0] != 4 and dic_reads[i].cigar[-1][1] >= 10:
			right_boundary.append(dic_reads[i].reference_end)
	if len(left_boundary) > 0:
		left_boundary_count = Counter(left_boundary)
		left_boundary_pos, left_boundary_freq = left_boundary_count.most_common(1)[0]
	if len(right_boundary) > 0:
		right_boundary_count = Counter(right_boundary)
		right_boundary_pos, right_boundary_freq = right_boundary_count.most_common(1)[0]
	# check each reads
	results = {}
	for i in dic_reads.keys():
		clipped = 0
		middle = 0
		non = 0
		if dic_reads[i].cigar[0][0] == 4 and dic_reads[i].cigar[-1][0] == 4 and dic_reads[i].cigar[0][1] >= 10 and dic_reads[i].cigar[-1][1] < 10:
			if left_boundary_pos - 5 <= dic_reads[i].reference_start and dic_reads[i].reference_start <= left_boundary_pos + 5:
				clipped += 1
		elif dic_reads[i].cigar[0][0] == 4 and dic_reads[i].cigar[-1][0] == 4 and dic_reads[i].cigar[0][1] < 10 and dic_reads[i].cigar[-1][1] >= 10:
			if right_boundary_pos - 5 <= dic_reads[i].reference_end and dic_reads[i].reference_end <= right_boundary_pos + 5:
				clipped += 1
		elif dic_reads[i].cigar[0][0] == 4 and dic_reads[i].cigar[-1][0] != 4 and dic_reads[i].cigar[0][1] >= 10:
			if left_boundary_pos - 5 <= dic_reads[i].reference_start and dic_reads[i].reference_start <= left_boundary_pos + 5:
				clipped += 1	
		elif dic_reads[i].cigar[-1][0] == 4 and dic_reads[i].cigar[0][0] != 4 and dic_reads[i].cigar[-1][1] >= 10:
			if right_boundary_pos - 5 <= dic_reads[i].reference_end and dic_reads[i].reference_end <= right_boundary_pos + 5:
				clipped += 1
		else:
			for cigartype, cigarlength in dic_reads[i].cigartuples:
				if cigartype == 1 and length_inserted_sequence-25 < cigarlength and cigarlength < length_inserted_sequence+25:
					middle += 1
		if clipped == 0 and middle == 0:
			non += 1
		results[i] = [clipped, middle, non]
	df_results = pd.DataFrame.from_dict(results, orient="index")
	if df_results.shape[0] == 0:
		vaf = -1
		return [0, 0, 0, -1]
	else:
		vaf = (df_results.loc[df_results[0]>0].shape[0] + df_results.loc[df_results[1]>0].shape[0]) / df_results.shape[0]
		return [df_results.loc[df_results[0]>0].shape[0], df_results.loc[df_results[1]>0].shape[0], df_results.shape[0], vaf]
#	return df_results


repeat_types = ["l1", "alu", "sva"]
celllines = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]

for repeat_type in repeat_types:
	for cellline in celllines:
		print(repeat_type, cellline)
		df = pd.read_csv("./results/%s_%s_tiers.txt"%(repeat_type, cellline), sep="\t")
		df_dip = pd.read_csv("./results/%s_%s_dip.txt"%(repeat_type, cellline), sep="\t")
		df_mini = pd.read_csv("./results/%s_%s_mini.txt"%(repeat_type, cellline), sep="\t")
		df_sni = pd.read_csv("./results/%s_%s_sni.txt"%(repeat_type, cellline), sep="\t")
		df_dip = df_dip.set_index("name_dip")
		df_mini = df_mini.set_index("name_mini")
		df_sni = df_sni.set_index("name_sni")
		df = df.loc[df["tier3"]=="tier3"]
		tmp = []
		for i in range(df.shape[0]):
			if df.iloc[i]["name_dip"] != ".":
				tmp.append(get_vaf(df.iloc[i]["name_dip"], df_dip, cellline))
			elif df.iloc[i]["name_sni"] != ".":
				tmp.append(get_vaf(df.iloc[i]["name_sni"], df_sni, cellline))
			else:
				tmp.append(get_vaf(df.iloc[i]["name_mini"], df_mini, cellline))
		df_final = pd.concat([df, pd.DataFrame(tmp, index=df.index, columns=["clipped", "middle", "total", "vaf"])], axis=1)
		df_final.sort_values(by="vaf")
		df_final.to_csv("%s_%s_tiers_vaf.txt"%(repeat_type, cellline), sep="\t", index=False)




