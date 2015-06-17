import os

data_path = "/backup/SNP/GRCh37_human"
index_path = data_path + "/bwa-index"
reads_path = data_path + "/chr20_mutate1/reads/4-1-2014"
results_path = data_path + "/chr20_mutate1/results/4-1-2014-BWA"

prog_path = "/home/nsvo/genome-tools/bwa-0.7.7"

for cvr in ['50x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = prog_path + "/bwa bwasw " + index_path + "/GRCh37_chr20.fasta " \
		+ reads_path + "/reads-100." + seq_err + "." + cvr + ".fq > " \
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.sam &"
		print cmd
		os.system(cmd)

