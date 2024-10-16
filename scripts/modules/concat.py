import pandas as pd




def concat_files(df_dip, df_mini, df_sni, repeat_type, cellline):
        df_dip_tmp = df_dip[["chrom", "pos", "name_dip"]]
        df_mini_tmp = df_mini[["chrom", "pos", "name_mini"]]
        df_sni_tmp = df_sni[["chrom", "pos", "name_sni"]]
        df_dip_tmp = df_dip_tmp.sort_values(by="pos")
        df_mini_tmp = df_mini_tmp.sort_values(by="pos")
        df_sni_tmp = df_sni_tmp.sort_values(by="pos")
        df_concat = pd.merge_asof(df_dip_tmp, df_sni_tmp, on='pos', by='chrom', tolerance=30, direction='nearest')
        df_sni_tmp['is_matched'] = df_sni_tmp['name_sni'].isin(df_concat['name_sni'])
        df_sni_tmp_unmatched= df_sni_tmp[df_sni_tmp['is_matched'] == False].drop(columns=['is_matched'])
        df_concat = pd.concat([df_concat, df_sni_tmp_unmatched], ignore_index=True)
        df_concat = df_concat.sort_values(by="pos")
        df_concat = pd.merge_asof(df_concat, df_mini_tmp, on='pos', by='chrom', tolerance=30, direction='nearest')
        df_mini_tmp['is_matched'] = df_mini_tmp['name_mini'].isin(df_concat['name_mini'])
        df_mini_tmp_unmatched= df_mini_tmp[df_mini_tmp['is_matched'] == False].drop(columns=['is_matched'])
        df_concat = pd.concat([df_concat, df_mini_tmp_unmatched], ignore_index=True)
        df_concat = df_concat.fillna(".")
        df_concat = df_concat.sort_values(by=["chrom", "pos"]).reset_index(drop=True)
        df_dip = df_dip.set_index("name_dip")
        df_mini = df_mini.set_index("name_mini")
        df_sni = df_sni.set_index("name_sni")
        info_dips = []
        info_minis = []
        info_snis = []
        for i in range(df_concat.shape[0]):
                if df_concat.iloc[i]["name_dip"] != ".":
                        info_dips.append(get_info(df_dip, df_concat.iloc[i]["name_dip"], repeat_type))
                else:
                        info_dips.append(".")
                if df_concat.iloc[i]["name_mini"] != ".":
                        info_minis.append(get_info(df_mini, df_concat.iloc[i]["name_mini"], repeat_type))
                else:
                        info_minis.append(".")
                if df_concat.iloc[i]["name_sni"] != ".":
                        info_snis.append(get_info(df_sni, df_concat.iloc[i]["name_sni"], repeat_type))
                else:
                        info_snis.append(".")
        df_concat["info_dip"] = info_dips
        df_concat["info_sni"] = info_snis
        df_concat["info_mini"] = info_minis
        df_concat = tiering(df_concat)
        df_concat = merge_information(df_concat, repeat_type, cellline)
        return df_concat



def merge_information(df_concat, repeat_type, cellline):
        total_blast = []
        total_tsd = []
        total_polyA = []
        total_direction = []
        for i in range(df_concat.shape[0]):
            if df_concat.iloc[i]["tier1"] == "tier1":
                blast, tsd, polyA, direction = extract_information(df_concat, i, "tier1", repeat_type, cellline)
                total_blast.append(blast)
                total_tsd.append(tsd)
                total_polyA.append(polyA)
                total_direction.append(direction)
            elif df_concat.iloc[i]["tier2"] == "tier2":
                blast, tsd, polyA, direction = extract_information(df_concat,i,  "tier2", repeat_type, cellline)
                total_blast.append(blast)
                total_tsd.append(tsd)
                total_polyA.append(polyA)
                total_direction.append(direction)
            elif df_concat.iloc[i]["tier3"] == "tier3":
                blast, tsd, polyA, direction = extract_information(df_concat, i, "tier3", repeat_type, cellline)
                total_blast.append(blast)
                total_tsd.append(tsd)
                total_polyA.append(polyA)
                total_direction.append(direction)
            else:
                total_blast.append(".")
                total_tsd.append(".")
                total_polyA.append(".")
                total_direction.append(".")
        df_concat["blast(subfamily;match_ratio;insert_start;insert_end;te_start;te_end)"] = total_blast
        df_concat["tsd"] = total_tsd
        df_concat["polyA"] = total_polyA
        df_concat["direction"] = total_direction
        return df_concat


def extract_blast_information(repeat_type, cellline, caller, name):
    df = pd.read_csv("../results/%s_%s_%s.txt"%(repeat_type, cellline, caller), sep="\t")
    blast = df.loc[df["name_%s"%caller] == name]["blast"].values[0]
    return blast


def extract_information(df_concat, i, tier, repeat_type, cellline):
    if df_concat.iloc[i]["info_dip"] != "." and df_concat.iloc[i]["info_dip"].split(";")[3] == tier:
        blast = extract_blast_information(repeat_type, cellline, "dip", df_concat.iloc[i]["name_dip"])
        tsd = df_concat.iloc[i]["info_dip"].split(";")[1]
        polyA = df_concat.iloc[i]["info_dip"].split(";")[2]
        if polyA == "polyA":
            direction = str(1)
        elif polyA == "polyT":
            direction = str(-1)
        else:
            direction = "."
        return blast, tsd, polyA, direction
    elif df_concat.iloc[i]["info_sni"] != "." and df_concat.iloc[i]["info_sni"].split(";")[3] == tier:
        blast = extract_blast_information(repeat_type, cellline, "sni", df_concat.iloc[i]["name_sni"])
        tsd = df_concat.iloc[i]["info_sni"].split(";")[1]
        polyA = df_concat.iloc[i]["info_sni"].split(";")[2]
        if polyA == "polyA":
            direction = str(1)
        elif polyA == "polyT":
            direction = str(-1)
        else:
            direction = "."
        return blast, tsd, polyA, direction
    elif df_concat.iloc[i]["info_mini"] != "." and df_concat.iloc[i]["info_mini"].split(";")[3] == tier:
        blast = extract_blast_information(repeat_type, cellline, "mini", df_concat.iloc[i]["name_mini"])
        tsd = df_concat.iloc[i]["info_mini"].split(";")[1]
        polyA = df_concat.iloc[i]["info_mini"].split(";")[2]
        if polyA == "polyA":
            direction = str(1)
        elif polyA == "polyT":
            direction = str(-1)
        else:
            direction = "."
        return blast, tsd, polyA, direction






def get_tiers(info, repeat_type):
        subfamily = info.split(";")[0]
        tsd = info.split(";")[1]
        polyA = info.split(";")[2]
        info_subfamily = 0
        info_tsd = 0
        info_polyA = 0
        if repeat_type == "l1":
                if subfamily == "L1HS" or subfamily == "L1PA2":
                        info_subfamily += 1
        elif repeat_type == "alu":
                if subfamily == "AluYa5" or subfamily == "AluYb8" or subfamily == "AluY":
                        info_subfamily += 1
        elif repeat_type == "sva":
                if subfamily == "SVA_F" or subfamily == "SVA_E":
                        info_subfamily += 1
        if tsd != ".":
                tsd_size = tsd.split("(")[0]
                if int(tsd_size) > 3 and int(tsd_size) < 31:
                        info_tsd += 1
        if polyA == "polyA" or polyA == "polyT":
                info_polyA += 1
        if info_subfamily > 0 and info_tsd > 0 and info_polyA > 0:
                return "tier1"
        elif info_subfamily > 0 and info_polyA > 0:
                return "tier2"
        elif info_subfamily > 0:
                return "tier3"
        else:
                return "."


def get_info(df, name, repeat_type):
        tmp = df.loc[name]
        if tmp["blast"] != ".;.;.;.;.;.":
                subfamily = tmp["blast"].split(";")[0]
        else:
                subfamily = "."
        if tmp["tsd"].split(";")[0] == "tsd":
                tsd = tmp["tsd"].split(";")[1]
        else:
                tsd = "."
        if tmp["polyA"].split(";")[0] == "polyA" or tmp["polyA"].split(";")[0] == "polyT":
                polyA = tmp["polyA"].split(";")[0]
        else:
                polyA = "."
        info = subfamily+";"+tsd+";"+polyA
        tiers = get_tiers(info, repeat_type)
        info = info+";"+tiers
        return info



def tiering(df):
        tier1 = []
        tier2 = []
        tier3 = []
        for i in range(df.shape[0]):
                tmp = []
                if df.iloc[i]["info_dip"] != ".":
                        if df.iloc[i]["info_dip"].split(";")[3] != ".":
                                tmp.append(int(df.iloc[i]["info_dip"].split(";")[3][-1]))
                if df.iloc[i]["info_mini"] != ".":
                        if df.iloc[i]["info_mini"].split(";")[3] != ".":
                                tmp.append(int(df.iloc[i]["info_mini"].split(";")[3][-1]))
                if df.iloc[i]["info_sni"] != ".":
                        if df.iloc[i]["info_sni"].split(";")[3] != ".":
                                tmp.append(int(df.iloc[i]["info_sni"].split(";")[3][-1]))
                if len(tmp) == 0:
                        tier1.append(".")
                        tier2.append(".")
                        tier3.append(".")
                else:
                        if min(tmp) == 1:
                                tier1.append("tier1")
                                tier2.append("tier2")
                                tier3.append("tier3")
                        elif min(tmp) == 2:
                                tier1.append(".")
                                tier2.append("tier2")
                                tier3.append("tier3")
                        elif min(tmp) == 3:
                                tier1.append(".")
                                tier2.append(".")
                                tier3.append("tier3")
        df["tier1"] = tier1
        df["tier2"] = tier2
        df["tier3"] = tier3
        return df

