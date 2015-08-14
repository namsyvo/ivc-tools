import os

data_path = "/backup/SNP/GRCh37_human"
reads_path = data_path + "/chr20_mutate1/reads/4-1-2014"

for cvr in ['50x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run convert_to_fq.go " + reads_path + "/reads-100." + seq_err + "." + cvr + ".txt > " \
		+ reads_path + "/reads-100." + seq_err + "." + cvr + ".fq &"
		print cmd
		os.system(cmd)

