"""
Get alignment info for fp (snp, indel)
Input: result folder, length of read/ref/mut
Ouput: alignment info (read-ref)
Usage: python get_diff_var_gatk_ivc.py cov_num var_type ivc_call_var_dir
"""
import sys
import os
import json

if len(sys.argv) != 5:
    print "Usage: python get_diff_var_gatk_ivc.py config_file cov_num var_type(0:tp_snp_unknown, 1:tp_snp_known, 2:tp_indel_unknown, 3:tp_indel_known) ivc_call_var_dir"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

data_path = data["DataPath"]["DataDir"]
result_dir = data["DataPath"]["ResultDir"]
ref_len = data["RefLen"]

cov_num = int(sys.argv[2])
var_type = int(sys.argv[3])
ivc_result_path = sys.argv[4]

read_lens = 100
read_nums = str(cov_num*ref_len/(2*read_lens))

'''
Find variants called by GATK but not by IVC
'''
ivc_fp_fn = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.20.0.0.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.20.0.0.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.20.0.0.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.20.0.0.0.txt"]

ivc_result_dn = os.path.join(data_path, result_dir, "ivc_0.70", ivc_result_path, "fpfntp_info")
var_pos = {}
fp_inf = open(os.path.join(ivc_result_dn, ivc_fp_fn[var_type]))
line = fp_inf.readline()
for line in fp_inf:
    tmp = line.strip().split()
    var_pos[tmp[0]] = True

gatk_fp_fn = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.20.0.txt", \
                  "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.20.0.txt", \
                  "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.20.0.txt", \
                  "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.20.0.txt"]

gatk_result_dn = os.path.join(data_path, result_dir, "gatk_hc_realign", "fpfntp_info")
diff_var_outf = open(os.path.join(gatk_result_dn, gatk_fp_fn[var_type] + "-diff-var-gatk-hc-not-" + ivc_result_path), "w")
fp_inf = open(os.path.join(gatk_result_dn, gatk_fp_fn[var_type]))
line = fp_inf.readline()
for line in fp_inf:
    tmp = line.strip().split()
    if tmp[0] not in var_pos:
        diff_var_outf.write(line)
diff_var_outf.close()

'''
Find variants called by IVC but not by GATK
'''
var_pos = {}
fp_inf = open(os.path.join(gatk_result_dn, gatk_fp_fn[var_type]))
line = fp_inf.readline()
for line in fp_inf:
    tmp = line.strip().split()
    var_pos[tmp[0]] = True

diff_var_outf = open(os.path.join(ivc_result_dn, ivc_fp_fn[var_type] + "-diff-var-ivc-not-gatk-hc"), "w")
diff_var_unique_outf = open(os.path.join(ivc_result_dn, ivc_fp_fn[var_type] + "-diff-var-unique_ivc-not-gatk-hc"), "w")
fp_inf = open(os.path.join(ivc_result_dn, ivc_fp_fn[var_type]))
line = fp_inf.readline()
ivc_var_pos = {}
prev_pos = ""
for line in fp_inf:
    tmp = line.strip().split()
    if tmp[0] not in var_pos:
        ivc_var_pos[tmp[0]] = True
        diff_var_outf.write(line)
        if tmp[0] != prev_pos:
            diff_var_unique_outf.write(line)
        prev_pos = tmp[0]
diff_var_outf.close()
diff_var_unique_outf.close()
