import os

data_path = "/home/SNP/GRCh37_human"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr20_mutate1/results/3-24-2014-BWA"
results_path = data_path + "/chr20_mutate1/results/GATK-time"

tool_path = "/home/nsvo/genome-tools"
script_path = "/home/nsvo/randalx/src/snpcaller/gatktools"

for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		cmd = "(time " + script_path + "/gatk-callsnp.sh " + ref_path + "/GRCh37_chr20.fasta " \
		+ sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted_RG.bam " \
		+ ref_path + "/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes-diffcontigname.vcf " \
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.vcf " + tool_path + ") 2>>gatk-hm-time.txt"
		os.system(cmd)
