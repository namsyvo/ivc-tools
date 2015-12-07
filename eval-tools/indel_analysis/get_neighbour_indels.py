'''
check_adjacent_indels
Usage: Usage: python check_adjacent_indels.py config_file prof_para input_file
'''

import os
import sys
import json
from datetime import datetime

if len(sys.argv) != 5:
    print "Usage: python check_adjacent_indels.py config_file prof_para input_file"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()
data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
dbsnp_fn = data["DataPath"]["dbsnpFile"]

para = sys.argv[2]
perf_meas = sys.argv[3]
input_file = sys.argv[4]

print "Getting ref genome and var prof info..."
var_prof = {}
chr_pos, chr_name = [], []
data = open(os.path.join(data_dir, index_dir, "index_0.70", "GRCh37.fasta.mgf.idx")).readlines()
for line in data:
    if line[0] == '>':
        info = line.strip().split()
        chr_name.append(info[0][1:])
        chr_pos.append(int(info[1]))
    else:
        break

var_prof = {}
for line in open(dbsnp_fn):
    if line[0] != '#':
        info = line.strip().split()
        offset_pos = -1
        for i in range(len(chr_pos)):
            if info[0] == chr_name[i]:
                offset_pos = chr_pos[i]
                break
        if offset_pos == -1:
            print "Missing chromosome", info[0]
        var_prof[offset_pos + int(info[1]) - 1] = info[3:5]

print "Getting true variants info..."
ref_path = os.path.join(data_dir, ref_dir)
known_var_file = os.path.join(ref_path, "known_var_" + para + ".txt")
unknown_var_file = os.path.join(ref_path, "unknown_var_" + para + ".txt")

true_known_snp, true_known_indel, true_unknown_snp, true_unknown_indel = {}, {}, {}, {}
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

print "Getting and evaluating called variants info..."
with open(input_file) as f:
    for line in f.readlines():
        value = line.strip().split()
        if len(value) >= 3 and value[1] == value[2]:
            continue
        try:
            var_pos = int(value[0])
            for p in range(30, 31):
	            k, u = 0, 0
	            for pos in range(var_pos - p, var_pos + p + 1):
	                if pos == var_pos:
	                    continue
	                if pos in true_known_indel:
	                    k += 1
	                elif pos in true_unknown_indel:
	                    u += 1
	            if k+u > 0:
	                print str(p) + "\t" + str(float(k)/(k + u)) + "\t" + perf_meas
	            else:
	            	print str(p) + "\t0\t" + perf_meas
        except ValueError:
            continue
