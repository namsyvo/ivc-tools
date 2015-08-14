import os

reads_path = "/home/SNP/GRCh37_human/chr20_mutate1/reads/3-24-2014"

for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run convert_to_fq.go " + reads_path + "/reads-100." + seq_err + "." + cvr + ".txt > " \
		+ reads_path + "/reads-100." + seq_err + "." + cvr + ".fq &"
		print cmd
		os.system(cmd)

