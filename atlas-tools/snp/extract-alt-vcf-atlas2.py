import os

data_path = "/backup/qmtran/Atlas_SNP2"
#result_path = data_path + "/human"
result_path = data_path
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['3x']:
	for seq_err in ['0.01']:
#for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "go run " + script_path + "/extract-alt-vcf.go " + result_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf.vcf > " \
		+ result_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf.txt &"
		os.system(cmd)
