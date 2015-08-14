import os

data_path = "/home/SNP/NC_016088.1_soybean"
result_path = data_path + "/chr1_mutate1/results/GATK-BOTH"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run " + script_path + "/extract-alt-vcf.go " + result_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.vcf > " \
		+ result_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.vcf.txt &"
		os.system(cmd)
