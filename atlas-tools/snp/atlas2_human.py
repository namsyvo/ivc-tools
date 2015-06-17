import os

data_path = "/backup/SNP/GRCh37_human"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr20_mutate1/results/3-24-2014-BWA"
results_path = "/backup/qmtran/Atlas_SNP2/human-time"

#script_path = "/home/SNP/Atlas2_v1.4.3/Atlas-SNP2"

for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.02']:

		cmd = "(time ruby Atlas-SNP2.rb -i " + sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted.bam " \
		+ "-r " + ref_path + "/GRCh37_chr20.fasta " + "-o "\
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf " \
		+ "-n GRCh37 --Illumina) 2>> atlas-snp2-hm-time.txt"
		print cmd
		os.system(cmd)
		'''
		cmd = "ruby Atlas-SNP2.rb -i " + sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted.bam " \
		+ "-r " + ref_path + "/GRCh37_chr20.fasta " + "-o "\
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf " \
		+ "-n GRCh37 -a " + ref_path + "/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes-diffcontigname.vcf -w --Illumina &"
		print cmd
		os.system(cmd)
		'''
