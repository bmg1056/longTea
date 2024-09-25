import os

def run_repeatmasker(repeat_type, cellline, caller, df):
    out = open("../tmp/repeatmasker/%s_%s_%s.fa"%(repeat_type, cellline, caller), 'w')
    for i in range(df.shape[0]):
        filename = df.iloc[i]["chrom"] + ";" + df.iloc[i]["pos"] + ";" + df.iloc[i]["caller"]
        out.write(">" + filename + "\n")
        out.write(df.iloc[i]["inserted_sequence"] + "\n")

    out.close()
    command = "RepeatMasker -species human ../tmp/repeatmasker/%s_%s_%s.fa"%(repeat_type, cellline, caller)
    os.system(command)




