
from Bio import SeqIO


fasta = list(SeqIO.parse("hg38_repeat.fa", "fasta"))

out = open("hg38_repeat.txt", 'w')
for i in fasta:
	out.write(i.id.upper() + "\t" + str(len(i.seq)) + "\n")

out.close()


