'''
Align reads with BWA-MEM
Usage: python bwa-align-exp.py config-dwgsim-chr1.txt 32 25 99cd33
'''
import os
import sys

config_file = sys.argv[1]
f=open(config_file)
prog_path = f.readline().strip()
data_path = f.readline().strip()
genome_fn = f.readline().strip()
read_fn = f.readline().strip()
f.close()

cpu_num = int(sys.argv[2])
cov_num = int(sys.argv[3])
rid = sys.argv[4]

ref_len = 249250621
read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [int(cov_num)]]

idx_path = data_path + "/indexes/bwa-index"
for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            read_path = data_path + "/reads/sim-reads/af_sid_mutant_dwgsim/alignment-analysis/"
            result_path = data_path + "/results/sim-reads/af_sid_mutant_dwgsim/bwa"
            if not os.path.exists(result_path):
                os.makedirs(result_path)
            #read_file_1 = read_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.read1.fastq." + rid
            #read_file_2 = read_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.read2.fastq." + rid
            read_file_1 = read_path + "/testread1." + rid
            read_file_2 = read_path + "/testread2." + rid
            result_file = result_path + "/bwa-analysis/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.sam." + rid
            cmd = prog_path + "/bwa mem -t " + str(cpu_num) + " " + idx_path + "/" + genome_fn + " " + \
                read_file_1 + " " + read_file_2 + " > " + result_file + " 2>" + result_file + ".log"
            print cmd
            os.system(cmd)
