data_path = "/home/genomes/SNP_calling"

read_path = data_path + "/reads/test-4M" + "/dag-gen-reads-checked-2-13-2014"
result_path = data_path + "/results/test-4M" + "/dag-gen-reads-checked-2-13-2014"

f = open(genome_path + "/eval.txt")
true_snp = {}
for line in f.readlines():
	if line.strip():
		value=line.strip().split()
		true_snp[value[0]]=value[1]

for cvr in ['5x']:
	for seq_err in ['0.01']:
		called_snp_file_path = result_path + "/called_snp-100." + seq_err + "." + cvr + ".txt"
		called_snp_file = open(called_snp_file_path)
		snp = {}
		for line in called_snp_file.readlines():
			if line.strip():
				value=line.strip().split()
				snp[value[0]]=value[1:]

		fp_snp_file_path = result_path + "/fp_snp-100." + seq_err + "." + cvr + ".txt"
		fp_snp_file = open(fp_snp_file_path, 'w')
		for key, value in snp.iteritems():
			if key in true_snp:
				if value[0] != true_snp[key]:
					fp_snp_file.write(str(value[0])+'\t'+str(value[1])+'\t'+str(value[2])\
					 +'\t'+str(value[3])+'\t'+str(true_snp[key])+'\t'+str(key)+'\n')
		fp_snp_file.close()
