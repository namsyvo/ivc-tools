import os

data_path = "/backup/SNP/GRCh37_human"
result_path = data_path + "/chr20_mutate1/results/4-1-2014-GATK"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['50x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run " + script_path + "/extract-alt-vcf.go " + result_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.vcf > " \
		+ result_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.vcf.txt &"
		os.system(cmd)
