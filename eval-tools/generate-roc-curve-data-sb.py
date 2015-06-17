import sys

data_path = "/backup/SNP/NC_016088.1_soybean"
index_path = data_path + "/index"
genome_path = data_path + "/chr1_mutate1/genome"
result_path = data_path + "/chr1_mutate1/results/RSC-ALL"

snp_file = index_path + "/SNPLocation.txt"
eval_file = genome_path + "/eval.txt"

true_snp = {}
true_indel = {}
with open(eval_file) as f:
	for line in f.readlines():
		if line.strip():
			value=line.strip().split()
			true_snp[int(value[0])] = value[1].strip()
			if len(value[1]) > 1:
				true_indel[int(value[0])] = value[1].strip()
confi_arr = [-50]
tmp = [n for n in range(50)]
confi_arr.extend(tmp)

for seq_err in ['0.02']:
	result_snp_file_path = result_path + "/prec-rec-snp-ISC-soybean-" + str(seq_err) + ".txt"
	result_snp_file = open(result_snp_file_path, "w")

	result_indel_file_path = result_path + "/prec-rec-indel-ISC-soybean-" + str(seq_err) + ".txt"
	result_indel_file = open(result_indel_file_path, "w")

	result_snp_file.write("Called-SNPs-seq_err-" + str(seq_err) + "\n")
	result_indel_file.write("Called-INDELs-seq_err-" + str(seq_err) + "\n")
	for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
		called_snp_file = result_path + "/called_snp-100." + seq_err + "." + cvr + ".txt"

		Prec_snp_arr = []
		Rec_snp_arr = []
		Prec_indel_arr = []
		Rec_indel_arr = []
		for confi in confi_arr:
			snp = {}
			with open(called_snp_file) as f:
				for line in f.readlines():
					if line.strip():
						value=line.strip().split()
						if len(value) == 5:
							if float(value[4]) >= (confi + 50)/100.0:
								snp[int(value[0])] = value[1].strip()
						else:
							if float(value[3]) >= (confi + 50)/100.0:
								snp[int(value[0])] = "."
			snp_num = len(true_snp)
			indel_num = len(true_indel)
			called_snp_num = 0
			true_called_snp_num = 0
			called_indel_num = 0
			true_called_indel_num = 0
			for key, value in snp.iteritems():
				if key in true_snp:
					called_snp_num += 1
					if len(value) > 1:
						called_indel_num += 1
					if value == true_snp[key]:
						true_called_snp_num += 1
						if len(value) > 1:
							true_called_indel_num += 1
			if called_snp_num != 0 and snp_num != 0:
				Prec_snp_arr.append(float(true_called_snp_num)/float(called_snp_num))
				Rec_snp_arr.append(float(true_called_snp_num)/float(snp_num - called_snp_num + true_called_snp_num))
			else:
				Prec_snp_arr.append(0.0)
				Rec_snp_arr.append(0.0)
			if called_indel_num != 0 and indel_num != 0:
				Prec_indel_arr.append(float(true_called_indel_num)/float(called_indel_num))
				Rec_indel_arr.append(float(true_called_indel_num)/float(indel_num - called_indel_num + true_called_indel_num))
			else:
				Prec_indel_arr.append(0.0)
				Rec_indel_arr.append(0.0)

		result_snp_file.write(cvr + "\n")
		result_indel_file.write(cvr + "\n")
		for confi in confi_arr:
			result_snp_file.write(str(Prec_snp_arr[confi]) + "\t")
			result_indel_file.write(str(Prec_indel_arr[confi]) + "\t")
		result_snp_file.write("\n")
		result_indel_file.write("\n")
		for confi in confi_arr:
			result_snp_file.write(str(Rec_snp_arr[confi]) + "\t")
			result_indel_file.write(str(Rec_indel_arr[confi]) + "\t")
		result_snp_file.write("\n")
		result_indel_file.write("\n")

	result_snp_file.close()
	result_indel_file.close()
