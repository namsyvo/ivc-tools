import os

data_path = "/home/SNP/GRCh37_human"
sam_path = data_path + "/chr20_mutate1/results/3-24-2014-BWA"

#for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
for cvr in ['1x']:
	for seq_err in ['0.01']:
		vcf_file = sam_path + "/reads-100." + seq_err + "." + cvr + ".bwasw.sam"
		f = open(vcf_file)
		lines = f.readlines()
		unaligned_num = 0
		prev_read_name = ""
		repeat_read_num = 0
		for line in lines:
			read_name = line.strip().split()[0]
			if read_name == prev_read_name:
				repeat_read_num += 1
				print read_name
			prev_read_name = line.strip().split()[0]
			aligned_flag = line.strip().split()[1]
			if aligned_flag == '4':
				unaligned_num += 1
		print unaligned_num
		print repeat_read_num
