import os

ref_path = "/home/SNP/GRCh37_human/genomes"
sam_path = "/home/SNP/GRCh37_human/chr20_mutate1/results/3-24-2014-BWA"
tool_path = "/home/nsvo/genome-tools"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

#for cvr in ['1x']:
#	for seq_err in ['0.01']:
for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = script_path + "/gatk-preprocess.sh " + ref_path + "/GRCh37_chr20 " \
		+ sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw " + tool_path + " &"
		os.system(cmd)
