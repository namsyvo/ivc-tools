'''
Evaluate variant call results for real data
'''

import os
import sys
import json
from datetime import datetime

if len(sys.argv) != 9:
    print "Usage: python eval_var_af_sid_mutant_diff_ref.py config_file confi_K_S confi_K_I confi_N_S confi_N_I result_dir_name"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

prog_ver = data["ProgVer"]
data_dir = data["DataPath"]["DataDir"]
dbsnp_fn = data["DataPath"]["dbsnpFile"]

confi_KS = float(sys.argv[2])
confi_KI = float(sys.argv[3])
confi_NS = float(sys.argv[4])
confi_NI = float(sys.argv[5])
result_dn = sys.argv[6]
cmp_tool_fn1 = sys.argv[7]
cmp_tool_fn2 = sys.argv[8]

var_prof = {}
var_prof_file = os.path.join(data_dir, "refs", dbsnp_fn)
KS, KI = 0, 0
with open(var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            var_prof[int(value[1]) - 1] = value[3:5]
            if len(value[3]) == 1 and len(value[4]) == 1:
                KS += 1
            else:
                KI += 1

cmp_tool_var1 = {}
with open(cmp_tool_fn1) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            if float(value[5]) >= 50:
                cmp_tool_var1[int(value[1]) - 1] = value[3:5]

cmp_tool_var2 = {}
with open(cmp_tool_fn2) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            if float(value[5]) >= 50:
                cmp_tool_var2[int(value[1]) - 1] = value[3:5]

result_path = os.path.join(data_dir, "results/real-reads", result_dn, "ivc")

var_call_file = os.path.join(result_path, "SRR352199.vcf")
var_call = {}
with open(var_call_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            if value[3] == value[4]:
                continue
            pos = int(value[1]) - 1
            if value[5] == "NaN":
                var_call[pos] = value[3:5]
            if pos in var_prof:
                if len(value[3]) == 1 and len(value[4]) == 1:
                    if float(value[5]) >= confi_KS and int(value[12]) > 1 and float(value[12])/float(value[13]) >= 0.25:
                        var_call[pos] = value[3:5]
                else:
                    if float(value[5]) >= confi_KI and int(value[12]) > 1 and float(value[12])/float(value[13]) >= 0.2:
                        var_call[pos] = value[3:5]
            else:
                if len(value[3]) == 1 and len(value[4]) == 1:
                    if float(value[5]) >= confi_NS: # and int(value[12]) > 1 and float(value[12])/float(value[13]) >= 0.25:
                        var_call[pos] = value[3:5]
                else:
                    if float(value[5]) >= confi_NI and int(value[12]) > 1 and float(value[12])/float(value[13]) >= 0.2:
                        var_call[pos] = value[3:5]

print "#called variants", len(var_call)

CS_KS, NC_KS, CS_KI, NC_KI, NS, NI = 0, 0, 0, 0, 0, 0
NS_CMP1, NI_CMP1, NS_CMP2, NI_CMP2, NS_CMP12, NI_CMP12 = 0, 0, 0, 0, 0, 0

for var_pos, var in var_call.iteritems():
    if var_pos in var_prof:
        if len(var[0]) == 1 and len(var[1]) == 1:
            if var == var_prof[var_pos]:
                CS_KS += 1
            else:
                #print var_pos, var, var_prof[var_pos]
                NC_KS += 1
        else:
            if var == var_prof[var_pos]:
                CS_KI += 1
            else:
                #print var_pos, var, var_prof[var_pos]
                NC_KI += 1
    else:
        if len(var[0]) == 1 and len(var[1]) == 1:
            NS += 1
            if var_pos in cmp_tool_var1 and var == cmp_tool_var1[var_pos]:
                NS_CMP1 += 1
            if var_pos in cmp_tool_var2 and var == cmp_tool_var2[var_pos]:
                NS_CMP2 += 1
            if var_pos in cmp_tool_var1 and var == cmp_tool_var1[var_pos] and var_pos in cmp_tool_var2 and var == cmp_tool_var2[var_pos]:
                NS_CMP12 += 1
        else:
            NI += 1
            if var_pos in cmp_tool_var1 and var == cmp_tool_var1[var_pos]:
                NI_CMP1 += 1
            if var_pos in cmp_tool_var2 and var == cmp_tool_var2[var_pos]:
                NI_CMP2 += 1
            if var_pos in cmp_tool_var1 and var == cmp_tool_var1[var_pos] and var_pos in cmp_tool_var2 and var == cmp_tool_var2[var_pos]:
                NI_CMP12 += 1

S = CS_KS + NC_KS + NS
I = CS_KI + NC_KI + NI
print "S, I, KS, KI", S, I, KS, KI
print "NS_CMP1, NI_CMP1, NS_CMP2, NI_CMP2, NS_CMP12, NI_CMP12", NS_CMP1, NI_CMP1, NS_CMP2, NI_CMP2, NS_CMP12, NI_CMP12

result_file_path = os.path.join(result_path, "SRR352199.prec_rec_time_mem.diff_pos.txt")
result_file = open(result_file_path, "w")
header = ["Alg", "Cov", "Qual", "CS_KS", "NC_KS", "CS_KI", "NC_KI", "NS", "NI", "S", "I", "V", \
              "CS_KS_ratio", "NC_KS_ratio", "CS_KI_ratio", "NC_KI_ratio", "NS_ratio", "NI_ratio", \
              "aln_num", "read_num", "aln_ratio", "memI", "timeI", "memC", "timeC", "memO", "timeO", "input_files", "input_paras", "prog_paras"]
result_file.write("\t".join(header) + "\n")
result_file.write(prog_ver + "\t" + str(75639683*2*100/3095677412.0) + "\t" + str(confi_KS) + "," + str(confi_KI) + "," + str(confi_NS) + "," + str(confi_NI) + "\t")
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t" % (CS_KS, NC_KS, CS_KI, NC_KI, NS, NI, S, I, S + I))
result_file.write("%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t" % (CS_KS/float(KS), NC_KS/float(KS), CS_KI/float(KI), NC_KI/float(KI), NS/float(S), NI/float(I)))

mem_time_file = os.path.join(result_path, "SRR352199.vcf.log")
with open(mem_time_file) as f:
    for line in f:
        tokens = line.strip().split("\t")
        if "Number of no-aligned reads" in tokens[0]:
            result_file.write(str(75639683-int(tokens[1])) + "\t75639683\t")
            result_file.write("%.5f\t" % ((75639683- int(tokens[1]))/float(75639683)))

with open(mem_time_file) as f:
    for line in f:
        tokens = line.strip().split("\t")
        if "Memstats after initializing the variant caller" in tokens[0]:
            result_file.write(tokens[3] + "\t")
        if "Time for initializing the variant caller" in tokens[0]:
            result_file.write(tokens[1] + "\t")
        if "Memstats after calling variants" in tokens[0]:
            result_file.write(tokens[3] + "\t")
        if "Time for calling variants" in tokens[0]:
            result_file.write(tokens[1] + "\t")
        if "Memstats after outputing variant calls" in tokens[0]:
            result_file.write(tokens[3] + "\t")
        if "Time for outputing variant calls" in tokens[0]:
            result_file.write(tokens[1] + "\t")
with open(mem_time_file) as f:
    for line in f:
        tokens = line.strip().split("\t")
        if "Input files" in tokens[0]:
            result_file.write(tokens[1] + "\t")
        if "Input paras" in tokens[0]:
            result_file.write(tokens[1] + "\t")
        if "Prog paras" in tokens[0]:
            result_file.write(tokens[1] + "\t")

result_file.close()
print "Check results at:", result_file_path
