'''
Call variants with IVC
Usage: python call_var_af_sid_mutant.py config_file coverage_num
'''
import os
import sys
import json
from datetime import datetime

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

prog_version = data["ProgVer"]
prog_path = data["ProgPath"]
data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
ref_fn = data["DataPath"]["GenomeFile"]
var_fn = data["DataPath"]["VarProfFile"]
read_dir = data["DataPath"]["ReadDir"]
index_dir = data["DataPath"]["IndexDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]
ref_len = data["RefLen"]

cov_num = sys.argv[2]
rid = ""
if len(sys.argv) == 4:
	rid = "." + sys.argv[3]
time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "_")
time_stamp = time_stamp.replace(":", "_")

ref_path = os.path.join(data_dir, ref_dir)
ref_file = os.path.join(ref_path, ref_fn)

ref_para = ['0.70', '0.75', '0.80', '0.85', '0.90', '0.95']
seq_errs = ['0.00015-0.0015']
read_lens = [100]
read_nums = []
if cov_num == "all":
    #read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 100]]
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25]]
else:
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [int(cov_num)]]

for para in ref_para[:1]:
    read_path = os.path.join(data_dir, read_dir)
    result_path = os.path.join(data_dir, result_dir, "ivc_" + para, prog_version + "_" + time_stamp)
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    var_prof_file = os.path.join(ref_path, var_fn + "_" + para + ".vcf")
    idx_dir = os.path.join(data_dir, index_dir, "index_" + para)
    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                    prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
                    if rid == "":
                            read_file_1 = os.path.join(read_path, prefix_fn + ".bwa.read1.fastq")
                            read_file_2 = os.path.join(read_path, prefix_fn + ".bwa.read2.fastq")
                    else:
                            read_file_1 = os.path.join(read_path, "diff_var_gatk_hc_analysis", prefix_fn + ".bwa.read1.fastq" + rid)
                            read_file_2 = os.path.join(read_path, "diff_var_gatk_hc_analysis", prefix_fn + ".bwa.read2.fastq" + rid)

                    var_call_file = os.path.join(result_path, prefix_fn + ".varcall.vcf")
                    time_mem_file = os.path.join(result_path, prefix_fn + ".varcall.log")
                    ivc_info_file = os.path.join(result_path, prefix_fn + ".varcall.info")

                    cmd = "(go run " + prog_path + " -R " + ref_file + " -V " + var_prof_file + " -I " + idx_dir + " -O " + var_call_file + \
                        " -1 " + read_file_1 + " -2 " + read_file_2 + " -debug true) 2>" + time_mem_file
                    print cmd
                    os.system(cmd)

cmd = "python eval_var_diff_ref_af_sid_mutant.py " + sys.argv[1] + " 20 7 1001 " + cov_num + " " + prog_version + "_" + time_stamp
print cmd
os.system(cmd)
