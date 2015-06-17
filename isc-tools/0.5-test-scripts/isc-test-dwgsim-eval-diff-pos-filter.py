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
genome_fn = data["DataPath"]["GenomeFile"]
snp_fn = data["DataPath"]["SNPProfFile"]
read_dir = data["DataPath"]["ReadDir"]
index_dir = data["DataPath"]["IndexDir"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]

confi = float(sys.argv[2])
cpu_num = sys.argv[3]
cov_num = sys.argv[4]
result_dn = sys.argv[5]

ref_path = os.path.join(data_dir, ref_dir)
read_path = os.path.join(data_dir, read_dir)

genome_file = os.path.join(ref_path, genome_fn)
snp_file = os.path.join(ref_path, snp_fn)

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

    true_snp_comp, true_indel_comp, true_snp_none, true_indel_none, all_var_prof, known_var_prof = {}, {}, {}, {}, {}, {}
    variant_comp_file = os.path.join(ref_path, "variant_comp_" + para + ".txt")
    variant_none_file = os.path.join(ref_path, "variant_none_" + para + ".txt")
    known_variant_prof_file = os.path.join(ref_path, "isc_snp_prof_" + para + ".vcf")
    all_variant_prof_file = os.path.join(ref_path, "snp_prof.vcf")

    with open(known_variant_prof_file) as f:
        for line in f.readlines():
            if line.strip() and line[0] != "#":
                value = line.strip().split()
                known_var_prof[int(value[1]) - 1] = value[3]

    with open(all_variant_prof_file) as f:
        for line in f.readlines():
            if line.strip() and line[0] != "#":
                value = line.strip().split()
                all_var_prof[int(value[1]) - 1] = value[3]

    KAKS, NS = 0, 0
    KAKID, NID = 0, 0

    with open(variant_comp_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                pos, val = int(value[0]), value[1]
                if len(val) == 1 and val != ".":
                    true_snp_comp[pos] = val
                    if val != all_var_prof[pos]:
                        KAKS += 1
                else:
                    true_indel_comp[pos] = val
                    if val != all_var_prof[pos]:
                        KAKID += 1

    with open(variant_none_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                pos, val = int(value[0]), value[1]
                if len(val) == 1 and val != ".":
                    true_snp_none[pos] = val
                    if val != all_var_prof[pos]:
                        NS += 1
                else:
                    true_indel_none[pos] = val
                    if val != all_var_prof[pos]:
                        NID += 1

    result_path = os.path.join(data_dir, result_dir, "isc_" + para, result_dn)
    result_file_path = result_path + "/" + read_fn + "_" + str(read_lens[0]) + "." + str(seq_errs[0]) + ".prec-rec-time-mem." + str(confi) + ".diff-pos.txt"
    result_file = open(result_file_path, "w")

    header = ["Alg", "cov", "qual", "TP-S", "FP-S", "FP-S-N", "TP-S-U", "FP-S-U", "TP-S-K", "FP-S-K", \
                "TP-I", "FP-I", "FP-I-N", "TP-I-U", "FP-I-U", "TP-I-K", "FP-I-K", \
                "P-S", "R-S", "P-S-U", "R-S-U", "P-S-K", "R-S-K", "P-I", "R-I", "P-I-U", "R-I-U", "P-I-K", "R-I-K", \
                "S", "S-U", "S-K", "CS", "CS-S", "I", "I-U", "I-K", "CI", "CI-I", \
                "run", "read", "proc", "max_snum", "max_ps_num", "na_num", "na_ratio", \
                "timeI", "memI", "timeC", "memC", "input_para", "para"]

    result_file.write("\t".join(header))
    result_file.write("\n")

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                for ms in max_snum[6:7]:
                    prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + "." + str(ms)
                    called_snp_file = os.path.join(result_path, prefix_fn + ".snpcall." + str(cpu_num) + ".vcf")
                    called_snp = {}
                    with open(called_snp_file) as f:
                        for line in f.readlines():
                            if line.strip():
                                value = line.strip().split()
                                key = int(value[0]) - 1
                                if float(value[2]) >= confi:
                                    if key in known_var_prof:
                                        if value[1] != known_var_prof[key]:
                                            if len(value[1]) == 2 and value[1][0] == value[1][1]:
                                                continue
                                            else:
                                                called_snp[key] = value[1]
                                    else:
                                        if len(value[1]) == 2 and value[1][0] == value[1][1]:
                                            continue
                                        else:
                                            called_snp[key] = value[1]

                    print "# called SNP", len(called_snp)

                    TP_KAKS, TP_NS = 0, 0
                    FP_KAKS, FP_NS = 0, 0
                    TP_KAKID, TP_NID = 0, 0
                    FP_KAKID, FP_NID = 0, 0
                    FP_S, FP_ID = 0, 0
                    same_none = 0
                    for key, value in called_snp.iteritems():
                        if key in true_snp_comp or key in true_indel_comp:
                            if key in true_snp_comp:
                                if value == true_snp_comp[key]:
                                    TP_KAKS += 1
                                else:
                                    FP_KAKS += 1
                            elif key in true_indel_comp:
                                if value == true_indel_comp[key]:
                                    TP_KAKID += 1
                                else:
                                    FP_KAKID += 1
                        elif key in true_snp_none or key in true_indel_none:
                            if key in true_snp_none:
                                if value == all_var_prof[key]:
                                    same_none += 1
                                if value == true_snp_none[key]:
                                    TP_NS += 1
                                else:
                                    FP_NS += 1
                            elif key in true_indel_none:
                                if value == true_indel_none[key]:
                                    TP_NID += 1
                                else:
                                    FP_NID += 1
                        else:
                            if len(value) == 1:
                                FP_S += 1
                            else:
                                FP_ID += 1
                    print "# known", TP_KAKS + FP_KAKS + TP_KAKID + FP_KAKID
                    print "# unknown", TP_NS + TP_NID, FP_NS + FP_NID, same_none
                    print "# other", FP_S + FP_ID

                    '''
                    nums_header = ["Alg", "cov", "qual", "TP-S", "FP-S", "FP-S-N", "TP-S-U", "FP-S-U", "TP-S-K", "FP-S-K", \
                                "TP-I", "FP-I", "FP-I-N", "TP-I-U", "FP-I-U", "TP-I-K", "FP-I-K"]
                    '''
                    result_file.write("\t".join([prog_version, "%.0f" % (2.0*int(rn)*int(rl)/ref_len), str(confi)]) + "\t")

                    #TP-S, FP-S, FP-S-N
                    result_file.write("%.5d\t" % (TP_KAKS + TP_NS))
                    result_file.write("%.5d\t" % (FP_KAKS + FP_NS + FP_S))
                    result_file.write(str(FP_S) + "\t")

                    #"TP-S-U", "FP-S-U", "TP-S-K", "FP-S-K"
                    result_file.write(str(TP_NS) + "\t" + str(FP_NS) + "\t")
                    result_file.write(str(TP_KAKS) + "\t" + str(FP_KAKS) + "\t")

                    #"TP-I", "FP-I", "FP-I-N"
                    result_file.write("%.5d\t" % (TP_KAKID + TP_NID))
                    result_file.write("%.5d\t" % (FP_KAKID + FP_NID + FP_ID))
                    result_file.write(str(FP_ID) + "\t")

                    #"TP-I-U", "FP-I-U", "TP-I-K", "FP-I-K"
                    result_file.write(str(TP_NID) + "\t" + str(FP_NID) + "\t")
                    result_file.write(str(TP_KAKID) + "\t" + str(FP_KAKID) + "\t")

                    '''
                    rate_header = ["P-S", "R-S", "P-S-U", "R-S-U", "P-S-K", "R-S-K", \
                                "P-I", "R-I", "P-I-U", "R-I-U", "P-I-K", "R-I-K"]
                    '''
                    #P-S, R-S
                    if TP_KAKS + FP_KAKS + TP_NS + FP_NS + FP_S != 0 and KAKS + NS != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKS + TP_NS)/float(TP_KAKS + TP_NS + FP_KAKS + FP_NS + FP_S)))
                        result_file.write("%.5f\t" % (float(TP_KAKS + TP_NS)/float(KAKS + NS)))
                    else:
                        result_file.write("\t\t")
                    #"P-S-U", "R-S-U"
                    if TP_NS + FP_NS != 0 and NS != 0:
                        result_file.write("%.5f\t" % (float(TP_NS)/float(TP_NS + FP_NS)))
                        result_file.write("%.5f\t" % (float(TP_NS)/float(NS)))
                    else:
                        result_file.write("\t\t")
                    #"P-S-K", "R-S-K"
                    if TP_KAKS + FP_KAKS != 0 and KAKS != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKS)/float(TP_KAKS + FP_KAKS)))
                        result_file.write("%.5f\t" % (float(TP_KAKS)/float(KAKS)))
                    else:
                        result_file.write("\t\t")
                    #"P-I", "R-I"
                    if TP_KAKID + FP_KAKID + TP_NID + FP_NID + FP_ID != 0 and KAKID + NID != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKID + TP_NID)/float(TP_KAKID + TP_NID + FP_KAKID + FP_NID + FP_ID)))
                        result_file.write("%.5f\t" % (float(TP_KAKID + TP_NID)/float(KAKID + NID)))
                    else:
                        result_file.write("\t\t")
                    #"P-I-U", "R-I-U"
                    if TP_NID + FP_NID != 0 and NID != 0:
                        result_file.write("%.5f\t" % (float(TP_NID)/float(TP_NID + FP_NID)))
                        result_file.write("%.5f\t" % (float(TP_NID)/float(NID)))
                    else:
                        result_file.write("\t\t")
                    #"P-I-K", "R-I-K"
                    if TP_KAKID + FP_KAKID != 0 and KAKID != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKID)/float(TP_KAKID + FP_KAKID)))
                        result_file.write("%.5f\t" % (float(TP_KAKID)/float(KAKID)))
                    else:
                        result_file.write("\t\t")

                    '''
                    para_header = ["S", "S-U", "S-K", "CS", "CS\S", "I", "I-U", "I-K", "CI", "CI\I", \
                                "run", "read", "proc", "max_snum", "max_ps_num", "na_num", "na_ratio", \
                                "timeI", "memI", "timeC", "memC", "input_para", "para"]
                    '''
                    #S, "S-U", "S-K", CS, CS-S
                    result_file.write(str(KAKS + NS) + "\t" + str(NS) + "\t" + str(KAKS) + "\t")
                    result_file.write("%.5d\t" % (FP_KAKS + FP_NS + FP_S + TP_KAKS + TP_NS))
                    result_file.write("%.5d\t" % ((FP_KAKS + FP_NS + FP_S + TP_KAKS + TP_NS) - (KAKS + NS)))

                    #"I", "I-U", "I-K", "CI", "CI-I"
                    result_file.write(str(KAKID + NID) + "\t" + str(NID) + "\t" + str(KAKID) + "\t")
                    result_file.write("%.5d\t" % (FP_KAKID + FP_NID + FP_ID + TP_KAKID + TP_NID))
                    result_file.write("%.5d\t" % ((FP_KAKID + FP_NID + FP_ID + TP_KAKID + TP_NID) - (KAKID + NID)))

                    #"run", "read", "proc", "max_snum", "max_ps_num"
                    result_file.write(result_dn + "\t" + prefix_fn + "\t" + cpu_num + "\t" + str(ms) + "\t1\t")
                    #"na_num", "na_ratio"
                    mem_time_file = os.path.join(result_path, prefix_fn + ".snpcall." + str(cpu_num) + ".log")
                    with open(mem_time_file) as f:
                        for line in f:
                            tokens = line.strip().split("\t")
                            if "# of no-aligned reads" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                                result_file.write(str((1-float(tokens[1]))/rn) + "\t")
                    #"timeI", "memI", "timeC", "memC"
                    with open(mem_time_file) as f:
                        for line in f:
                            tokens = line.strip().split("\t")
                            if "time for initializing SNP caller" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                            if "memstats after initializing SNP caller" in tokens[0]:
                                result_file.write(str(float(tokens[3])/10**9) + "\t")
                            if "time for calling SNPs" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                            if "memstats after calling SNPs" in tokens[0]:
                                result_file.write(str(float(tokens[3])/10**9) + "\t")
                    #"input_para", "para"
                    with open(mem_time_file) as f:
                        for line in f:
                            tokens = line.strip().split("\t")
                            if "Input parameters" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                            if "Parameters" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                    result_file.write("\n")

result_file.close()
