import os

data_path = "/backup/SNP/GRCh37_human"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr20_mutate1/results/4-1-2014-BWA"
tool_path = "/home/nsvo/genome-tools"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['50x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = script_path + "/gatk-preprocess.sh " + ref_path + "/GRCh37_chr20 " \
		+ sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw " + tool_path + " &"
		os.system(cmd)
