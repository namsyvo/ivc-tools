'''
Evaluate variant call results for real data
'''

import os
import sys
import json
from datetime import datetime

if len(sys.argv) != 6:
    print "Usage: python eval_var_af_sid_mutant_diff_ref.py config_file confi_S confi_I result_dir_name"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()
prog_ver = data["ProgVer"]
data_dir = data["DataPath"]["DataDir"]
dbsnp_fn = data["DataPath"]["dbsnpFile"]

confi_S = float(sys.argv[2])
confi_I = float(sys.argv[3])
result_dn = sys.argv[4]
cmp_tool_fn = sys.argv[5]

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

cmp_tool_var = {}
with open(cmp_tool_fn) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            if float(value[5]) >= 50:
                cmp_tool_var[int(value[1]) - 1] = value[3:5]

result_path = os.path.join(data_dir, "results/real-reads", result_dn, "gatk_hc_realign")

var_call_file = os.path.join(result_path, "SRR352199.bwa_sorted_RG_realign_Noknown.vcf")
var_call = {}
with open(var_call_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            if value[3] == value[4]:
                continue
            pos = int(value[1]) - 1
            if len(value[3]) == 1 and len(value[4]) == 1:
                if float(value[5]) >= confi_S:
                    var_call[pos] = value[3:5]
            else:
                if float(value[5]) >= confi_I:
                    var_call[pos] = value[3:5]

print "#called variants", len(var_call)

CS_KS, NC_KS, CS_KI, NC_KI, NS, NI = 0, 0, 0, 0, 0, 0
NS_CMP, NI_CMP = 0, 0
for var_pos, var in var_call.iteritems():
    if var_pos in var_prof:
        if len(var[0]) == 1 and len(var[1]) == 1:
            if var == var_prof[var_pos]:
                CS_KS += 1
            else:
                NC_KS += 1
        else:
            if var == var_prof[var_pos]:
                CS_KI += 1
            else:
                NC_KI += 1
    else:
        if len(var[0]) == 1 and len(var[1]) == 1:
            NS += 1
            if var_pos in cmp_tool_var and var == cmp_tool_var[var_pos]:
                NS_CMP += 1
        else:
            NI += 1
            if var_pos in cmp_tool_var and var == cmp_tool_var[var_pos]:
                NI_CMP += 1

S = CS_KS + NC_KS + NS
I = CS_KI + NC_KI + NI

result_file_path = os.path.join(result_path, "SRR352199.bwa_sorted_RG_realign_Noknown.vcf.prec_rec_time_mem.diff_pos.txt")
result_file = open(result_file_path, "w")

header = ["Alg", "Cov", "Qual", "CS_KS", "NC_KS", "CS_KI", "NC_KI", "NS", "NI", "NS_CMP", "NI_CMP", "S", "I", "V", \
              "CS_KS_ratio", "NC_KS_ratio", "CS_KI_ratio", "NC_KI_ratio", "NS_ratio", "NI_ratio", "time"]

result_file.write("\t".join(header) + "\n")
result_file.write(prog_ver + "\t" + str(75639683*2*100/3095677412.0) + "\t" + str(confi_S) + "," + str(confi_I) + "\t")
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t" % (CS_KS, NC_KS, CS_KI, NC_KI, NS, NI, NS_CMP, NI_CMP, S, I, S+I))
result_file.write("%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t" % (CS_KS/float(KS), NC_KS/float(KS), CS_KI/float(KI), NC_KI/float(KI), NS/float(S), NI/float(I)))

mem_time_file = os.path.join(result_path, "SRR352199.bwa_sorted_RG_realign_Noknown.vcf.log")
with open(mem_time_file) as f:
    for line in f:
        if "Total runtime" in line:
            result_file.write(line)

result_file.close()
print "Check results at:", result_file_path
