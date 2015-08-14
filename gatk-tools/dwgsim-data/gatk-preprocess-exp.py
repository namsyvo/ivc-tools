'''
Preprocess SAM files with GATK
Usage: python gatk-preprocess-exp.py config-dwgsim-chr1.txt
'''
import os
import sys

config_file = sys.argv[1]
f=open(config_file)
prog_path = f.readline().strip()
data_path = f.readline().strip()
genome_fn = f.readline().strip()
dbsnp_fn = f.readline().strip()
read_fn = f.readline().strip()
f.close()

ref_len = 249250621
read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [50, 100]]
#read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [2, 3, 4, 6, 7, 8, 9, 10, 15, 20, 25, 30]]

ref_path = data_path + "/refs"
tool_path = "/home/nsvo/genome-tools"

for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            sam_path = data_path + "/results/sim-reads/af_sid_mutant_dwgsim/bwa"
            cmd = prog_path + "/gatk-preprocess.sh " + ref_path + "/GRCh37_chr1 " \
                + sam_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa " + tool_path + " &"
            os.system(cmd)
