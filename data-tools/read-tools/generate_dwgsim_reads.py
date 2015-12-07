'''
Generate reads with DWGSIM
Usage: python dwgsim-genreads-exp.py
'''
import os

ref_len = 3095677412 #Whole Genome
#ref_len = 249250621 #chr1
#ref_len = 243199373 #chr2
#ref_len = 492449994 #chr1_2
#ref_len = 1062541960 #chr1_5
#ref_len = 1815907890 #chr1_10
#ref_len = 2881033286 #chr_all

read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [40, 50, 100]]
#read_nums = [100000]
for rl in read_lens:
	for err in seq_errs:
		for rn in read_nums:
			ref_path = "/data/nsvo/test_data/GRCh37/refs/af_sid_mutant/mutant_genome.fasta"
			result_path = "/data/nsvo/test_data/GRCh37/reads/sim_reads/af_sid_mutant_dwgsim"
			if not os.path.exists(result_path):
				os.makedirs(result_path)
			cmd = "/home/SNP/genome-tools/DWGSIM/dwgsim -N " + str(rn) + " -1 " + str(rl) + " -2 " + str(rl) \
			+ " -e " + str(err) + " -E " + str(err) + " -r 0 " + ref_path + " " + result_path \
			+ "/dwgsim_reads_" + str(rl) + "." + str(err) + "." + str(rn) + " &"
			print cmd
			os.system(cmd)
