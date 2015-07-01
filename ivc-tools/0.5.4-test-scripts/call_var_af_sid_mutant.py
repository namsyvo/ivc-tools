'''
Call variant from FASTQ file
Usage: python call_var_af_sid_mutant.py config_file cpu_num coverage_num
    E.g.: python call_var_af_sid_mutant.py config-chr1-test.json 32 5
'''
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
var_fn = data["DataPath"]["VarProfFile"]
read_dir = data["DataPath"]["ReadDir"]
index_dir = data["DataPath"]["IndexDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]

cpu_num = sys.argv[2]
cov_num = sys.argv[3]
rid = ""
if len(sys.argv) == 5:
	rid = "." + sys.argv[4]
time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "-")

ref_path = os.path.join(data_dir, ref_dir)
genome_file = os.path.join(ref_path, genome_fn)

ref_len = 249250621
ref_para = ['0.70', '0.75', '0.80', '0.85', '0.90', '0.95', '0.96', '0.97', '0.98', '0.99']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
max_snum = [2**i for i in range(3, 14)]
read_nums = []
if cov_num == "all":
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25]]
else:
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [int(cov_num)]]

for para in ref_para[0:1]:
    read_path = os.path.join(data_dir, read_dir)
    result_path = os.path.join(data_dir, result_dir, "ivc_" + para, prog_version + "-" + time_stamp)
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    var_file = os.path.join(ref_path, var_fn + "_" + para + ".vcf")
    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                for ms in max_snum[6:7]:
	                prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
	                #read_file_1 = os.path.join(read_path, prefix_fn + ".bwa.read1.fastq" + rid)
	                #read_file_2 = os.path.join(read_path, prefix_fn + ".bwa.read2.fastq" + rid)
	                read_file_1 = os.path.join(read_path, "alignment-analysis", prefix_fn + ".bwa.read1.fastq" + rid)
	                read_file_2 = os.path.join(read_path, "alignment-analysis", prefix_fn + ".bwa.read2.fastq" + rid)

	                mem_time_file = os.path.join(result_path, prefix_fn + "." + str(ms) + ".varcall." + str(cpu_num) + ".log")
	                info_file = os.path.join(result_path, prefix_fn + "." + str(ms) + ".varcall." + str(cpu_num) + ".info")
	                cpu_prof_file = os.path.join(result_path, prefix_fn + "." + str(ms) + ".varcall." + str(cpu_num) + ".pprof")
	                called_var_file = os.path.join(result_path, prefix_fn + "." + str(ms) + ".varcall." + str(cpu_num) + ".vcf")
	                cmd = "(go run " + prog_path + \
	                    " -g " + genome_file + " -s " + var_file + " -i " + os.path.join(data_dir, index_dir, "index_" + para) + \
	                    " -1 " + read_file_1 + " -2 " + read_file_2 + " -o " + called_var_file  + " -c " + cpu_prof_file + \
	                    " -w " + cpu_num + " -t " + cpu_num + " -n " + str(ms) + " -h 30 -k 128) 2>" + mem_time_file + " 1>" + info_file
	                print cmd
	                os.system(cmd)

#cmd = "python eval_var_diff_ref_af_sid_mutant.py " + config_file + " 1.14 " + cpu_num + " " + cov_num + " " + prog_version + "-" + time_stamp
#print cmd
#os.system(cmd)

