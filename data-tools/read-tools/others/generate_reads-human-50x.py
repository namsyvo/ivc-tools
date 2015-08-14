import os

data_path = "/backup/SNP/GRCh37_human"
genome_path = data_path + "/chr20_mutate1/genome"
reads_path = data_path + "/chr20_mutate1/reads/4-1-2014"

for cvr in ['100']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run generate_reads.go -c " + cvr + " -e " + seq_err + " -s " + genome_path + \
			"/mutate.txt > " + reads_path + "/reads-100." + seq_err + "." + cvr + "x.txt &"
		print cmd
		os.system(cmd)
