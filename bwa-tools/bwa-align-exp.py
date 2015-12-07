'''
Align reads with BWA-MEM
Usage: python bwa-align-exp.py config-dwgsim-chr1.txt 32
'''
import os
import sys
import json

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

prog_path = data["ProgPath"]
data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
ref_fn = data["DataPath"]["RefFile"]
read_dir = data["DataPath"]["ReadDir"]
index_dir = data["DataPath"]["IndexDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]
ref_len = data["RefLen"]

cpu_num = int(sys.argv[2])

read_lens = [100]
seq_errs = ['0.00015-0.0015']
#read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 100]]
#read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 5, 10, 15, 20, 25]]
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1]]

for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            read_path = data_dir + "/" + read_dir
            result_path = data_dir + "/" + result_dir + "/bwa"
            if not os.path.exists(result_path):
                os.makedirs(result_path)
            read_file_1 = read_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.read1.fastq"
            read_file_2 = read_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.read2.fastq"
            result_file = result_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.sam"
            cmd = "/usr/bin/time -v " + prog_path + "/bwa mem -t " + str(cpu_num) + " " + data_dir + "/" + index_dir + "/" + ref_fn + " " + \
                read_file_1 + " " + read_file_2 + " > " + result_file + " 2>" + result_file + ".log"
            print cmd
            os.system(cmd)
