import os

genome_path = "/home/SNP/GRCh37_human/chr20_mutate1/genome"
read_path = "/home/SNP/GRCh37_human/chr20_mutate1/reads/3-24-2014"

for cvr in ['1', '2', '3', '5', '7', '10']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run generate_reads.go -c " + cvr + " -e " + seq_err + " -s " + genome_path + \
			"/mutate.txt > " + read_path + "/reads-100." + seq_err + "." + cvr + "x.txt &"
		print cmd
		os.system(cmd)
