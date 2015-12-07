'''
Preprocess SAM files with GATK
Usage: python gatk-preprocess-exp.py config_file cov_num
'''
import os
import sys
import json

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

prog_path = data["ProgPath"]
script_path = data["ScriptPath"]
data_dir = data["DataPath"]["DataDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]
ref_len = data["RefLen"]

cov_num = sys.argv[2]

read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = []
if cov_num == "all":
    #read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 100]]
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [15, 20, 25]]
else:
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [int(cov_num)]]

ref_file = os.path.join(data_dir, "refs", "GRCh37")
for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            sam_path = os.path.join(data_dir, result_dir, "bwa")
            cmd = "/usr/bin/time -v " + script_path + "/gatk-preprocess.sh " + ref_file + " " \
                + sam_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa " + prog_path + \
                " 2>" + sam_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.preprocess.log"
            os.system(cmd)
