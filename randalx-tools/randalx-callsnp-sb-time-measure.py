import os

data_path = "/home/SNP/NC_016088.1_soybean"
ref_path = data_path + "/index"
read_path = data_path + "/chr1_mutate1/reads/all-reads"
result_path = data_path + "/chr1_mutate1/results/RSC-ALL"

prog_path = "/home/nsvo/randalx/src/randalx-align/randalx.go"

genome_file = ref_path + "/genomestar.txt"
snp_file = ref_path + "/SNPLocation.txt"
index_file = ref_path + "/genomestar.txt.index"
rev_index_file = ref_path + "/genomestar_rev.txt.index"

for cvr in ['1x', '2x', '3x', '5x', '7x', '10x']:
	for seq_err in ['0.01', '0.02', '0.04']:
		query_file = read_path + "/reads-100." + seq_err + "." + cvr + ".txt"
		called_snp_file = result_path + "/called_snp-100." + seq_err + "." + cvr + ".txt"

		print "Running with coverage ", cvr, " sequencing error ", seq_err
		cmd = "(time go run " + prog_path + " -g " + genome_file + " -s " + snp_file + \
			" -i " + index_file + " -r " + rev_index_file + " -q " + query_file + \
			" -c " + called_snp_file + " -e " + seq_err + " -l 100" + \
			" > " + result_path + "/prog_stat-100." + seq_err + "." + cvr + ".txt) 2>>rsc-sb-time.txt"

		print cmd
		os.system(cmd)
