import os
import sys
import json
from datetime import datetime

config_file = sys.argv[1]
f = open(config_file)
data = json.load(f)
f.close()

prog_version = data["ProgVer"]
prog_path = data["ProgPath"]
data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
genome_fn = data["DataPath"]["GenomeFile"]
snp_fn = data["DataPath"]["SNPProfFile"]
read_dir = data["DataPath"]["ReadDir"]
index_dir = data["DataPath"]["IndexDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]

cpu_num = sys.argv[2]

time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "-")

ref_path = os.path.join(data_dir, ref_dir)
genome_file = os.path.join(ref_path, genome_fn)
snp_file = os.path.join(ref_path, snp_fn)

ref_len = 249250621
ref_para = ['0.0000', '0.0825', '0.1650', '0.2475', '0.3300']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [5]]
max_snum = [2**i for i in range(3, 14)]

for para in ref_para[0:1]:
    read_path = os.path.join(data_dir, read_dir + "/mutate-" + para + "-dwgsim")
    result_path = os.path.join(data_dir, result_dir + "/mutate-" + para + "-dwgsim", "isc", "debug", prog_version + "-" + time_stamp)
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                for ms in max_snum[9:10]:
	                prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
	                read_file_1 = os.path.join(read_path, prefix_fn + ".bwa.read1.fastq")
	                read_file_2 = os.path.join(read_path, prefix_fn + ".bwa.read2.fastq")

	                mem_time_file = os.path.join(result_path, prefix_fn + "." + str(ms) + ".snpcall." + str(cpu_num) + ".log")
	                called_snp_file = os.path.join(result_path, prefix_fn + "." + str(ms) + ".snpcall." + str(cpu_num) + ".vcf")
	                cmd = "(go run " + prog_path + \
	                    " -g " + genome_file + " -s " + snp_file + " -i " + os.path.join(data_dir, index_dir) + \
	                    " -1 " + read_file_1 + " -2 " + read_file_2 + " -o " + called_snp_file + \
	                    " -w " + cpu_num + " -t " + cpu_num + " -n " + str(ms) + ") 2>" + mem_time_file
	                '''
	                cmd = "go run " + prog_path + \
	                    " -g " + genome_file + " -s " + snp_file + " -i " + os.path.join(data_dir, index_dir) + \
	                    " -1 " + read_file_1 + " -2 " + read_file_2 + " -o " + called_snp_file + \
	                    " -w " + cpu_num + " -t " + cpu_num
	                '''
	                print cmd
	                os.system(cmd)
	                print cmd

cmd = "python isc-test-dwgsim-eval-chr1.py " + config_file + " 24.83 " + cpu_num + " " + prog_version + "-" + time_stamp
os.system(cmd)
