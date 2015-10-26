'''
check_adjacent_indels
Usage: Usage: python check_adjacent_indels.py config_file prof_para input_file
'''

import os
import sys
import json
from datetime import datetime

if len(sys.argv) != 4:
    print "Usage: python check_adjacent_indels.py config_file prof_para input_file"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()
data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
dbsnp_fn = data["DataPath"]["dbsnpFile"]

para = sys.argv[2]
input_file = sys.argv[3]

var_prof = {}
var_prof_file = os.path.join(data_dir, "refs", dbsnp_fn)

with open(var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            var_prof[int(value[1]) - 1] = value[3:5]

ref_path = os.path.join(data_dir, ref_dir)
true_known_snp, true_known_indel, true_unknown_snp, true_unknown_indel = {}, {}, {}, {}
known_var_file = os.path.join(ref_path, "known_var_" + para + ".txt")
unknown_var_file = os.path.join(ref_path, "unknown_var_" + para + ".txt")

with open(known_var_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            pos, known_var = int(value[0]), value[1:]
            if len(var_prof[pos][0]) == 1 and len(var_prof[pos][1]) == 1:
                true_known_snp[pos] = known_var
            else:
                true_known_indel[pos] = known_var

with open(unknown_var_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            pos, unknown_var = int(value[0]), value[1:]
            if len(var_prof[pos][0]) == 1 and len(var_prof[pos][1]) == 1:
                true_unknown_snp[pos] = unknown_var
            else:
                true_unknown_indel[pos] = unknown_var

K_w, N_w_k, N_W_u = {}, {}, {}
var_num = 0
with open(input_file) as f:
    for line in f.readlines():
        value = line.strip().split()
        if len(value) >= 3 and value[1] == value[2]:
            continue
        try:
            var_pos = int(value[0])
            for p in range(1, 31):
                k, u = 0, 0
                for pos in range(var_pos - p, var_pos + p + 1):
                    if pos == var_pos:
                        continue
                    if pos in true_known_indel:
                        N_w_k[p] = N_w_k.get(p, 0) + 1
                        k += 1
                    elif pos in true_unknown_indel:
                        N_W_u[p] = N_W_u.get(p, 0) + 1
                        u += 1
                if k+u == 0:
                    K_w[p] = K_w.get(p, 0) + 0
                else:
                    K_w[p] = K_w.get(p, 0) + float(k)/(k + u)
            var_num += 1
        except ValueError:
            print "Error", ValueError, "at", line

print "\t#var\t" + str(var_num)
print "p\tn_i_k\tn_i_u\tn_i_sum\tk_ratio\tn_i_norm\tn_i_k_ratio"

for p in range(1, 31):
    try:
        N_w_sum = N_w_k[p] + N_W_u[p]
        print str(p) + "\t" + str(N_w_k[p]) + "\t" + str(N_W_u[p]) + "\t" + str(N_w_sum) + "\t" \
            + str(float(N_w_k[p])/N_w_sum) + "\t" + str(float(N_w_sum)/var_num) + "\t" + str(float(K_w[p])/var_num)
    except KeyError:
        print "Error", KeyError, "at", p
