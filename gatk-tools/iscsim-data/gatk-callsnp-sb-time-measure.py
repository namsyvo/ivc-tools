import os

data_path = "/home/SNP/NC_016088.1_soybean"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr1_mutate1/results/BWA-ALL-SAMS"
results_path = data_path + "/chr1_mutate1/results/GATK-time"

tool_path = "/home/nsvo/genome-tools"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "(time " + script_path + "/gatk-callsnp.sh " + ref_path + "/chr1.fasta " \
		+ sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted_RG.bam " \
		+ ref_path + "/vcf_chr_1-diffcontigname.vcf " \
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.vcf " + tool_path + ") 2>>gatk-sb-time.txt"
		os.system(cmd)
