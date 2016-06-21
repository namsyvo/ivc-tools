import os

prog_path = "/home/nsvo/randalx/src/randalx-align/randalx.go"
data_path = "/home/genomes/SNP_calling"

genome_path = data_path + "/genomes/Af293chr1"
read_path = data_path + "/reads/Af293chr1" + "/2-15-2014"
result_path = data_path + "/results/Af293chr1" + "/2-15-2014"

genome_file = genome_path + "/genomestar.txt"
snp_file = genome_path + "/SNPLocation.txt"
index_file = genome_path + "/genomestar.txt.index"
rev_index_file = genome_path + "/genomestar_rev.txt.index"
eval_file = genome_path + "/eval.txt"

#for cvr in ['5x', '10x', '20x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
#for cvr in ['5x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
for cvr in ['3x']:
	for seq_err in ['0.01']:
#for cvr in ['3x', '5x', '7x', '10x']:
#	for seq_err in ['0.01', '0.02', '0.04']:
		query_file = read_path + "/reads-100." + seq_err + "." + cvr + ".txt"
		#query_file = read_path + "/xaa"
		called_snp_file = result_path + "/test_called_snp-100." + seq_err + "." + cvr + ".txt"

		#cmd = "rm " + " ".join([path + "/called_snp_fp_check.txt", path + "/stat_fp_check.txt", \
		#	path + "/exact_nf.txt", path + "/extended_nf.txt", path + "/aligned.txt", path + "/matched_nf.txt"])
		#os.system(cmd)
		cmd = "go run " + prog_path + " -g " + genome_file + " -s " + snp_file + \
			" -i " + index_file + " -r " + rev_index_file + " -q " + query_file + \
			" -c " + called_snp_file + " -e " + seq_err + " -l 100" + \
			" > " + result_path + "/test_prog_stat-100." + seq_err + "." + cvr + ".txt &"

		'''
		cmd = "go run " + prog_path + " -g " + genome_file + " -s " + snp_file + \
			" -i " + index_file + " -r " + rev_index_file + " -q " + query_file + \
			" -c " + called_snp_file + " -e " + seq_err + " -l 100"
		'''

		print cmd
		os.system(cmd)
