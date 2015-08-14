import os

data_path = "/backup/SNP/NC_016088.1_soybean"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr1_mutate1/results/BWA-ALL-SAMS"
tool_path = "/home/nsvo/genome-tools"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['3x']:
	for seq_err in ['0.01']:
		cmd = script_path + "/gatk-preprocess.sh " + ref_path + "/chr1 " \
		+ sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw " + tool_path + " &"
		os.system(cmd)
