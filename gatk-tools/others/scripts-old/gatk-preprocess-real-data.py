import os

ref_path = "/backup/SNP/GRCh37_human/genomes"
sam_path = "/data/nsvo/Human/SRS000204/results"
tool_path = "/data/nsvo/genome-tools"
script_path = "/data/nsvo/gatk-scripts"

cmd = script_path + "/gatk-preprocess-real-data-1.sh " + ref_path + "/GRCh37_chr20 " \
	+ sam_path + "/SRR352199-chr1 " + tool_path
os.system(cmd)
