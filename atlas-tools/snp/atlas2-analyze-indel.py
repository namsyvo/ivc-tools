data_path = "/backup/SNP/GRCh37_human"
#data_path = "/backup/SNP/NC_016088.1_soybean"
index_path = data_path + "/index"
genome_path = data_path + "/chr20_mutate1/genome"
#genome_path = data_path + "/chr1_mutate1/genome"
result_path = "/backup/qmtran/Atlas_Indel2/human"
#result_path = "/backup/qmtran/Atlas_SNP2/soybean"

snp_file = index_path + "/SNPLocation.txt"
eval_file = genome_path + "/eval.txt"

true_snp = {}
true_indel = {}
with open(eval_file) as f:
	for line in f.readlines():
		if line.strip():
			value=line.strip().split()
			true_snp[int(value[0])]=value[1]
			if len(value[1]) > 1:
				true_indel[int(value[0])]=value[1]

called_indel_fn = result_path + "/called_indel-length-atlas2-human.txt"
called_indel_file = open(called_indel_fn, "w")
true_called_indel_fn = result_path + "/true_called_indel-length-atlas2-human.txt"
true_called_indel_file = open(true_called_indel_fn, "w")

for cvr in ['5x']:
	for seq_err in ['0.02']:
#for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
		called_snp_file = result_path + "/reads-100." + seq_err + "." + cvr + ".atlas2.vcf.txt"
		snp = {}
		with open(called_snp_file) as f:
			for line in f.readlines():
				if line.strip():
					value=line.strip().split()
					snp[int(value[0])]=value[1]
		for key, value in snp.iteritems():
			if key - 1 in true_snp:
				if len(value) > 1:
					called_indel_file.write(str(len(value)) + "\n")
				if value == true_snp[key - 1]:
					if len(value) > 1:
						true_called_indel_file.write(str(len(value)) + "\n")
called_indel_file.close()
true_called_indel_file.close()
