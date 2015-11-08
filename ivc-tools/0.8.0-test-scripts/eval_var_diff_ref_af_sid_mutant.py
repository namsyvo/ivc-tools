'''
Evaluate variant call results
Usage: python eval_var_af_sid_mutant_diff_ref.py config_file confi_K confi_U cov_num result_dir_name
    E.g.: python eval_var_af_sid_mutant_diff_ref.py config_chr1_test.json 3.0 20.0 5 IVC_0.5.3_aff_gap_aln_2015_05_20_15:09:36.226789
'''

import os
import sys
import json
from datetime import datetime

if len(sys.argv) != 11:
    print "Usage: python eval_var_af_sid_mutant_diff_ref.py config_file confi_K_S confi_K_I confi_U_S confi_U_I confi_1_S confi_1_I prof_para cov_num result_dir_name"
    exit(0)

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]
dbsnp_fn = data["DataPath"]["dbsnpFile"]
ref_len = data["RefLen"]

confi_K_S = float(sys.argv[2])
confi_K_I = float(sys.argv[3])
confi_U_S = float(sys.argv[4])
confi_U_I = float(sys.argv[5])
confi_1_S = float(sys.argv[6])
confi_1_I = float(sys.argv[7])
para = sys.argv[8]
cov_num = sys.argv[9]
result_dn = sys.argv[10]

seq_errs = ['0.00015-0.0015']
read_lens = [100]
read_nums = []
if cov_num == "all":
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 100]]
else:
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [int(cov_num)]]

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

KAKS, NS = 0, 0
KAKID, NID = 0, 0
with open(known_var_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            pos, known_var = int(value[0]), value[1:]
            if len(var_prof[pos][0]) == 1 and len(var_prof[pos][1]) == 1:
                true_known_snp[pos] = known_var
                if known_var[0] != known_var[1]:
                    KAKS += 1
            else:
                true_known_indel[pos] = known_var
                if known_var[0] != known_var[1]:
                    KAKID += 1

with open(unknown_var_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            pos, unknown_var = int(value[0]), value[1:]
            if len(var_prof[pos][0]) == 1 and len(var_prof[pos][1]) == 1:
                true_unknown_snp[pos] = unknown_var
                if unknown_var[0] != unknown_var[1]:
                    NS += 1
            else:
                true_unknown_indel[pos] = unknown_var
                if unknown_var[0] != unknown_var[1]:
                    NID += 1

result_path = os.path.join(data_dir, result_dir, "ivc_" + para, result_dn)
result_file_path = result_path + "/" + read_fn + "_" + str(read_lens[0]) + "." + str(seq_errs[0]) + \
    ".prec_rec_time_mem." + str(confi_K_S) + "." + str(confi_K_I) + "." + str(confi_U_S) + "." + str(confi_U_I) + "." + cov_num + ".diff_pos.txt"

result_file = open(result_file_path, "w")

header = ["Alg", "Cov", "Qual", "TP_S", "FP_S", "FP_S_N", "TP_S_U", "FP_S_U", "TP_S_K", "FP_S_K", \
            "TP_I", "FP_I", "FP_I_N", "TP_I_U", "FP_I_U", "TP_I_K", "FP_I_K", \
            "P_S", "R_S", "P_S_U", "R_S_U", "P_S_K", "R_S_K", "P_I", "R_I", "P_I_U", "R_I_U", "P_I_K", "R_I_K", \
            "S", "S_U", "S_K", "CS", "CS_S", "I", "I_U", "I_K", "CI", "CI_I", \
            "na_num", "na_ratio", "memI", "timeI", "memC", "timeC", "memO", "timeO", "input_files", "input_paras", "prog_paras"]

result_file.write("\t".join(header))
result_file.write("\n")

for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
            var_call_file = os.path.join(result_path, prefix_fn + ".varcall.vcf")
            var_call = {}
            cn = rn/(ref_len/(2*read_lens[0]))
            with open(var_call_file) as f:
                for line in f.readlines():
                    if line.strip() and line[0] != '#':
                        value = line.strip().split()
                        if value[3] == value[4]:
                            continue
                        pos = int(value[1]) - 1
                        if value[5] == "NaN":
                            var_call[pos] = value[3:5]
                        if pos in true_known_snp:
                            if float(value[5]) >= confi_K_S:
                                var_call[pos] = value[3:5]
                        elif pos in true_known_indel:
                            if float(value[5]) >= confi_K_I:
                                var_call[pos] = value[3:5]
                        elif len(value[3]) == 1 and len(value[4]) == 1:
                            if float(value[5]) >= confi_U_S:
                                if cn <= 5:
                                    if float(value[12]) > 1.0:
                                        var_call[pos] = value[3:5]
                                    elif float(value[5]) >= confi_1_S and float(value[13]) < 2:
                                        var_call[pos] = value[3:5]
                                else:
                                    if float(value[12])/cn > 0.25:
                                        var_call[pos] = value[3:5]
                        else:
                            if float(value[5]) >= confi_U_I:
                                if cn <= 5:
                                    if float(value[12]) > 1.0:
                                        var_call[pos] = value[3:5]
                                    elif float(value[5]) >= confi_1_I and float(value[13]) < 2:
                                        var_call[pos] = value[3:5]
                                else:
                                    if float(value[12])/cn > 0.20:
                                        var_call[pos] = value[3:5]
            print "#called variants", len(var_call)

            TP_KAKS, TP_NS = 0, 0
            FP_KAKS, FP_NS = 0, 0
            TP_KAKID, TP_NID = 0, 0
            FP_KAKID, FP_NID = 0, 0
            FP_S, FP_ID = 0, 0
            for var_pos, var in var_call.iteritems():
                if var_pos in true_known_snp or var_pos in true_known_indel:
                    if var_pos in true_known_snp:
                        if var == true_known_snp[var_pos]:
                            TP_KAKS += 1
                        else:
                            FP_KAKS += 1
                    elif var_pos in true_known_indel:
                        if var == true_known_indel[var_pos]:
                            TP_KAKID += 1
                        else:
                            FP_KAKID += 1
                elif var_pos in true_unknown_snp or var_pos in true_unknown_indel:
                    if var_pos in true_unknown_snp:
                        if var == true_unknown_snp[var_pos]:
                            TP_NS += 1
                        else:
                            FP_NS += 1
                    elif var_pos in true_unknown_indel:
                        if var == true_unknown_indel[var_pos]:
                            TP_NID += 1
                        else:
                            FP_NID += 1
                else:
                    if len(var[0]) == 1 and len(var[1]) == 1:
                        FP_S += 1
                    else:
                        FP_ID += 1

            print "# known", TP_KAKS + FP_KAKS + TP_KAKID + FP_KAKID
            print "# unknown", TP_NS + TP_NID, FP_NS + FP_NID
            print "# none", FP_S + FP_ID

            '''
            nums_header = ["run", "cov", "qual", "TP_S", "FP_S", "FP_S_N", "TP_S_U", "FP_S_U", "TP_S_K", "FP_S_K", \
                        "TP_I", "FP_I", "FP_I_N", "TP_I_U", "FP_I_U", "TP_I_K", "FP_I_K"]
            '''
            result_file.write("\t".join([result_dn, "%.0f" % (2.0*int(rn)*int(rl)/ref_len), str(confi_K_S) + ", " + \
                str(confi_K_I) + ", " + str(confi_U_S) + ", " + str(confi_U_I) + ", " + str(confi_1_S) + ", " + str(confi_1_I)]) + "\t")

            #TP_S, FP_S, FP_S_N
            result_file.write("%.5d\t" % (TP_KAKS + TP_NS))
            result_file.write("%.5d\t" % (FP_KAKS + FP_NS + FP_S))
            result_file.write(str(FP_S) + "\t")

            #"TP_S_U", "FP_S_U", "TP_S_K", "FP_S_K"
            result_file.write(str(TP_NS) + "\t" + str(FP_NS) + "\t")
            result_file.write(str(TP_KAKS) + "\t" + str(FP_KAKS) + "\t")

            #"TP_I", "FP_I", "FP_I_N"
            result_file.write("%.5d\t" % (TP_KAKID + TP_NID))
            result_file.write("%.5d\t" % (FP_KAKID + FP_NID + FP_ID))
            result_file.write(str(FP_ID) + "\t")

            #"TP_I_U", "FP_I_U", "TP_I_K", "FP_I_K"
            result_file.write(str(TP_NID) + "\t" + str(FP_NID) + "\t")
            result_file.write(str(TP_KAKID) + "\t" + str(FP_KAKID) + "\t")

            '''
            rate_header = ["P_S", "R_S", "P_S_U", "R_S_U", "P_S_K", "R_S_K", \
                        "P_I", "R_I", "P_I_U", "R_I_U", "P_I_K", "R_I_K"]
            '''
            #P_S, R_S
            if TP_KAKS + FP_KAKS + TP_NS + FP_NS + FP_S != 0 and KAKS + NS != 0:
                result_file.write("%.5f\t" % (float(TP_KAKS + TP_NS)/float(TP_KAKS + TP_NS + FP_KAKS + FP_NS + FP_S)))
                result_file.write("%.5f\t" % (float(TP_KAKS + TP_NS)/float(KAKS + NS)))
            else:
                result_file.write("\t\t")
            #"P_S_U", "R_S_U"
            if TP_NS + FP_NS != 0 and NS != 0:
                result_file.write("%.5f\t" % (float(TP_NS)/float(TP_NS + FP_NS + FP_S)))
                result_file.write("%.5f\t" % (float(TP_NS)/float(NS)))
            else:
                result_file.write("\t\t")
            #"P_S_K", "R_S_K"
            if TP_KAKS + FP_KAKS != 0 and KAKS != 0:
                result_file.write("%.5f\t" % (float(TP_KAKS)/float(TP_KAKS + FP_KAKS)))
                result_file.write("%.5f\t" % (float(TP_KAKS)/float(KAKS)))
            else:
                result_file.write("\t\t")
            #"P_I", "R_I"
            if TP_KAKID + FP_KAKID + TP_NID + FP_NID + FP_ID != 0 and KAKID + NID != 0:
                result_file.write("%.5f\t" % (float(TP_KAKID + TP_NID)/float(TP_KAKID + TP_NID + FP_KAKID + FP_NID + FP_ID)))
                result_file.write("%.5f\t" % (float(TP_KAKID + TP_NID)/float(KAKID + NID)))
            else:
                result_file.write("\t\t")
            #"P_I_U", "R_I_U"
            if TP_NID + FP_NID != 0 and NID != 0:
                result_file.write("%.5f\t" % (float(TP_NID)/float(TP_NID + FP_NID + FP_ID)))
                result_file.write("%.5f\t" % (float(TP_NID)/float(NID)))
            else:
                result_file.write("\t\t")
            #"P_I_K", "R_I_K"
            if TP_KAKID + FP_KAKID != 0 and KAKID != 0:
                result_file.write("%.5f\t" % (float(TP_KAKID)/float(TP_KAKID + FP_KAKID)))
                result_file.write("%.5f\t" % (float(TP_KAKID)/float(KAKID)))
            else:
                result_file.write("\t\t")

            '''
            para_header = ["S", "S_U", "S_K", "CS", "CS\S", "I", "I_U", "I_K", "CI", "CI\I", \
                        "na_num", "na_ratio", "memI", "timeI", "memC", "timeC", "memO", "timeO", \
                        "input_files", "input_paras", "prog_paras"]
            '''
            #S, "S_U", "S_K", CS, CS_S
            result_file.write(str(KAKS + NS) + "\t" + str(NS) + "\t" + str(KAKS) + "\t")
            result_file.write("%.5d\t" % (FP_KAKS + FP_NS + FP_S + TP_KAKS + TP_NS))
            result_file.write("%.5d\t" % ((FP_KAKS + FP_NS + FP_S + TP_KAKS + TP_NS) - (KAKS + NS)))

            #"I", "I_U", "I_K", "CI", "CI_I"
            result_file.write(str(KAKID + NID) + "\t" + str(NID) + "\t" + str(KAKID) + "\t")
            result_file.write("%.5d\t" % (FP_KAKID + FP_NID + FP_ID + TP_KAKID + TP_NID))
            result_file.write("%.5d\t" % ((FP_KAKID + FP_NID + FP_ID + TP_KAKID + TP_NID) - (KAKID + NID)))

            #"na_num", "na_ratio"
            mem_time_file = os.path.join(result_path, prefix_fn + ".varcall.log")
            with open(mem_time_file) as f:
                for line in f:
                    tokens = line.strip().split("\t")
                    if "Number of no-aligned reads" in tokens[0]:
                        result_file.write(tokens[1] + "\t")
                        result_file.write(str((1-float(tokens[1]))/rn) + "\t")
            #"timeI", "memI", "timeC", "memC", "timeO", "memO"
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
            #"input_para", "para"
            with open(mem_time_file) as f:
                for line in f:
                    tokens = line.strip().split("\t")
                    if "Input files" in tokens[0]:
                        result_file.write(tokens[1] + "\t")
                    if "Input paras" in tokens[0]:
                        result_file.write(tokens[1] + "\t")
                    if "Prog paras" in tokens[0]:
                        result_file.write(tokens[1] + "\t")
            result_file.write("\n")
result_file.close()
print "Check results at:", result_file_path
