data_path = "/home/genomes/SNP_calling"

genome_path = data_path + "/genomes/Af293chr1"
result_path = data_path + "/results/Af293chr1" + "/2-15-2014"

snp_file = genome_path + "/SNPLocation.txt"
eval_file = genome_path + "/eval.txt"

true_snp = {}
with open(eval_file) as f:
	for line in f.readlines():
		if line.strip():
			value=line.strip().split()
			true_snp[int(value[0])]=value[1]

result_file_path = result_path + "/test_prec_rec_all.txt"
result_file = open(result_file_path, "w")

stat_file_path = result_path + "/test_stat_all.txt"
stat_file = open(stat_file_path, "w")

#for cvr in ['5x', '10x', '20x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
#for cvr in ['5x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
for cvr in ['3x']:
	for seq_err in ['0.01']:

#for cvr in ['3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
		called_snp_file = result_path + "/test_called_snp-100." + seq_err + "." + cvr + ".txt"
		snp = {}
		with open(called_snp_file) as f:
			for line in f.readlines():
				if line.strip():
					value=line.strip().split()
					snp[int(value[0])]=value[1]

		snp_num = len(true_snp)
		called_snp_num = 0
		true_called_snp_num = 0
		for key, value in snp.iteritems():
			if key in true_snp:
				called_snp_num += 1
				if value == true_snp[key]:
					true_called_snp_num += 1

		result_file.write(cvr + "\t" + seq_err + "\t100\n")
		result_file.write("snp_num\t" + str(snp_num) + "\n")
		result_file.write("called_snp_num\t" + str(called_snp_num) + "\n")
		result_file.write("true_called_snp_num\t" + str(true_called_snp_num) + "\n")

		if called_snp_num != 0 and snp_num != 0:
			result_file.write("Precision\t" + str(float(true_called_snp_num)/float(called_snp_num)) + "\n")
			result_file.write("Recall\t" + str(float(true_called_snp_num)/float(snp_num - called_snp_num + true_called_snp_num)) + "\n")
		else:
			result_file.write("Precision\t\n")
			result_file.write("Recall\t\n")
		result_file.write("\n")

		stat_path = result_path + "/test_prog_stat-100." + seq_err + "." + cvr + ".txt"
		with open(stat_path) as f:
			lines = f.readlines()

		stat_file.write(cvr + "\t" + seq_err + "\t100\n")
		stat_file.write(lines[1])
		stat_file.write(lines[2])
		stat_file.write(lines[3])
		stat_file.write(lines[5])
		stat_file.write(lines[6])

result_file.close()
stat_file.close()
