import os

data_path = "/backup/SNP/NC_016088.1_soybean"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr1_mutate1/results/BWA-ALL-SAMS"
results_path = "/backup/qmtran/Atlas_SNP2/"

#script_path = "/home/SNP/Atlas2_v1.4.3/Atlas-SNP2"

for cvr in ['3x']:
	for seq_err in ['0.01']:
#for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
		'''		
		cmd = "(time ruby Atlas-SNP2.rb -i " + sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted.bam " \
		+ "-r " + ref_path + "/chr1.fasta " + "-o "\
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf " \
		+ "-n NC_016088.1 -a " + ref_path + "/vcf_chr_1-diffcontigname.vcf --Illumina) 2>> atlas-snp2-sb-time.txt"
		print cmd
		os.system(cmd)
		'''
		cmd = "ruby Atlas-SNP2.rb -i " + sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted.bam " \
		+ "-r " + ref_path + "/chr1.fasta " + "-o "	+ results_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf " \
		+ "-n NC_016088.1 --Illumina"
		print cmd
		os.system(cmd)

