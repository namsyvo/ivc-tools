"""
Get alignment info for fp (snp, indel)
Input: result folder, length of read/ref/mut
Ouput: alignment info (read-ref)
Usage: python get_diff_var_gatk_ivc.py cov_num var_type ivc_call_var_dir

0: tp_snp_unknown
1: tp_snp_known
2: tp_indel_unknown
3: tp_indel_known
4: fp_snp_known
5: fp_snp_unknown
6: fp_snp_none
7: fp_indel_known
8: fp_indel_unknown
9: fp_indel_none

"""
import sys
import os
import json

if len(sys.argv) != 6:
    print "Usage: python get_diff_var_gatk_ivc.py config_file cov_num var_type call_var_dir1 call_var_dir2"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

data_path = data["DataPath"]["DataDir"]
result_dir = data["DataPath"]["ResultDir"]
ref_len = data["RefLen"]

cov_num = int(sys.argv[2])
var_type = int(sys.argv[3])
result_path1 = sys.argv[4]
result_path2 = sys.argv[5]

read_lens = 100
read_nums = str(cov_num*ref_len/(2*read_lens))

'''
Find variants called by GATK but not by IVC
'''
file_name = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".fp_snp_known.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".fp_snp_unknown.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".fp_snp_none.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".fp_indel_known.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".fp_indel_unknown.20.0.7.0.txt", \
                 "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".fp_indel_none.20.0.7.0.txt"]

result_dn1 = os.path.join(data_path, result_dir, "ivc_0.70", result_path1, "fpfntp_info")
var_pos1 = {}
fp_inf = open(os.path.join(result_dn1, file_name[var_type]))
line = fp_inf.readline()
for line in fp_inf:
    tmp = line.strip().split()
    var_pos1[tmp[0]] = True
fp_inf.close()

result_dn2 = os.path.join(data_path, result_dir, "ivc_0.70", result_path2, "fpfntp_info")
var_pos2 = {}
fp_inf = open(os.path.join(result_dn2, file_name[var_type]))
line = fp_inf.readline()
for line in fp_inf:
    tmp = line.strip().split()
    var_pos2[tmp[0]] = True
fp_inf.close()

'''
Find variants called by 1 but not 2
'''
diff_var_outf = open(os.path.join(result_dn1, file_name[var_type] + "-diff-var-not-" + result_path2), "w")
diff_var_unique_outf = open(os.path.join(result_dn1, file_name[var_type] + "-diff-var-unique-not-" + result_path2), "w")
fp_inf = open(os.path.join(result_dn1, file_name[var_type]))
line = fp_inf.readline()
prev_pos = ""
for line in fp_inf:
    tmp = line.strip().split()
    if tmp[0] not in var_pos2:
        diff_var_outf.write(line)
        if tmp[0] != prev_pos:
            diff_var_unique_outf.write(line)
        prev_pos = tmp[0]
diff_var_outf.close()
diff_var_unique_outf.close()

'''
Find variants called by 2 but not 1
'''
diff_var_outf = open(os.path.join(result_dn2, file_name[var_type] + "-diff-var-not-" + result_path1), "w")
diff_var_unique_outf = open(os.path.join(result_dn2, file_name[var_type] + "-diff-var-unique-not-" + result_path1), "w")
fp_inf = open(os.path.join(result_dn2, file_name[var_type]))
line = fp_inf.readline()
ivc_var_pos = {}
prev_pos = ""
for line in fp_inf:
    tmp = line.strip().split()
    if tmp[0] not in var_pos1:
        ivc_var_pos[tmp[0]] = True
        diff_var_outf.write(line)
        if tmp[0] != prev_pos:
            diff_var_unique_outf.write(line)
        prev_pos = tmp[0]
diff_var_outf.close()
diff_var_unique_outf.close()
