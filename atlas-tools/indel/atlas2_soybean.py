import os

data_path = "/backup/SNP/NC_016088.1_soybean"
ref_path = data_path + "/genomes"
sam_path = data_path + "/chr1_mutate1/results/BWA-ALL-SAMS"
results_path = "/backup/qmtran/Atlas_Indel2/soybean-test"

#script_path = "/backup/SNP/Atlas2_v1.4.3/Atlas-SNP2"

for cvr in ['3x']:
	for seq_err in ['0.01']:
#for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
		'''
		cmd = "(time ruby Atlas-Indel2.rb -b " + sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted.bam " \
		+ "-r " + ref_path + "/chr1.fasta " + "-o "\
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".atlas2-indel.vcf " \
		+ "-I -a " + ref_path + "/vcf_chr_1-diffcontigname.vcf) 2>> atlas-indel2-sb-time.txt"
		print cmd
		os.system(cmd)
		'''
		cmd = "ruby Atlas-Indel2.rb -b " + sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw_sorted.bam " \
		+ "-r " + ref_path + "/chr1.fasta " + "-o "\
		+ results_path + "/reads-100." + seq_err + "." + cvr + ".atlas2-indel.vcf " \
		+ "-I"
		print cmd
		os.system(cmd)
