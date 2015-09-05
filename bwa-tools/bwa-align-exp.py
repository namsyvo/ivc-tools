'''
Align reads with BWA-MEM
Usage: python bwa-align-exp.py config-dwgsim-chr1.txt 32
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

#ref_len = 249250621 #chr1
ref_len = 243199373 #chr2
read_lens = [100]
seq_errs = ['0.00015-0.0015']
#read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 100]]
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [10]]

idx_path = data_path + "/indexes/bwa-index"
for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            read_path = data_path + "/reads/sim-reads/af_sid_mutant_dwgsim"
            result_path = data_path + "/results/sim-reads/af_sid_mutant_dwgsim/bwa"
            if not os.path.exists(result_path):
                os.makedirs(result_path)
            read_file_1 = read_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.read1.fastq"
            read_file_2 = read_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.read2.fastq"
            result_file = result_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.sam"
            cmd = prog_path + "/bwa mem -t " + str(cpu_num) + " " + idx_path + "/" + genome_fn + " " + \
                read_file_1 + " " + read_file_2 + " > " + result_file + " 2>" + result_file + ".log"
            print cmd
            os.system(cmd)
