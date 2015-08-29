import os
import sys
import json
import re
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
result_dn = sys.argv[4]

ref_path = os.path.join(data_dir, ref_dir)
read_path = os.path.join(data_dir, read_dir)

genome_file = os.path.join(ref_path, genome_fn)
snp_file = os.path.join(ref_path, snp_fn)

ref_len = 249250621
ref_para = ['0.70', '0.75', '0.80', '0.85', '0.90', '0.95', '0.96', '0.97', '0.98', '0.99']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
#read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25]]
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [5]]
max_snum = [2**i for i in range(3, 14)]

for para in ref_para[0:1]:

    true_snp_comp, true_indel_comp, true_snp_none, true_indel_none = {}, {}, {}, {}

    variant_comp_file = os.path.join(ref_path, "alt_mutant", "variant_comp_" + para + ".txt")
    variant_none_file = os.path.join(ref_path, "alt_mutant", "variant_none_" + para + ".txt")

    with open(variant_comp_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                if len(value[1]) == 1 and value[1] != ".":
                    true_snp_comp[int(value[0])] = value[1]
                else:
                    true_indel_comp[int(value[0])] = value[1]

    with open(variant_none_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                if len(value[1]) == 1 and value[1] != ".":
                    true_snp_none[int(value[0])] = value[1]
                else:
                    true_indel_none[int(value[0])] = value[1]

    KAKS, NS = len(true_snp_comp), len(true_snp_none)
    KAKID, NID = len(true_indel_comp), len(true_indel_none)

    result_path = os.path.join(data_dir, result_dir, "alt_mutant_dwgsim", "isc_" + para, result_dn)
    result_file_path = result_path + "/" + read_fn + "_" + str(read_lens[0]) + "." + str(seq_errs[0]) + ".prec-rec-time-mem.dis-filter." + str(confi) + ".txt"
    result_file = open(result_file_path, "w")

    header = ["Alg", "cov", "qual", \
                "S", "P-S", "R-S", "FP-S@Other", "P-S@None", "R-S@None", "P-S@Comp", "R-S@Comp", \
                "I", "P-I", "R-I", "FP-I@Other", "P-I@None", "R-I@None", "P-I@Comp", "R-I@Comp", \
                "S@None", "TP-S@None", "FP-S@None", "S@Comp", "TP-S@Comp", "FP-S@Comp", \
                "I@None", "TP-I@None", "FP-I@None", "I@Comp", "TP-I@Comp", "FP-I@Comp", \
                "run", "read", "proc", "max_snum", "max_ps_num", "a_num", "na_num", \
                "timeI", "memI", "timeC", "memC", "input_para", "para"]

    result_file.write("\t".join(header))
    result_file.write("\n")

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                for ms in max_snum[6:7]:
                    prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + "." + str(ms)
                    called_snp_file = os.path.join(result_path, prefix_fn + ".snpcall." + str(cpu_num) + ".vcf")
                    snp = {}
                    with open(called_snp_file) as f:
                        for line in f.readlines():
                            if line.strip():
                                value = line.strip().split()
                                if float(value[3]) >= confi:
                                    snp[int(value[0]) - 1] = value[1:]

                    sorted_snp = sorted(snp.iteritems(), key=lambda x: float(x[1][2]), reverse=True)
                    print "# called snps after filter1", len(sorted_snp)

                    called_snp = {}
                    expected_snp_num = 2895847
                    snp_num = 0
                    for key, value in sorted_snp:
                        if key in true_snp_comp or key in true_indel_comp:
                            called_snp[key] = value
                            if len(value[0]) == 1:
                                snp_num += 1

                    for key, value in sorted_snp:
                        aln_dis = int(value[4])
                        chr_dis = int(value[5])
                        if key in true_snp_comp or key in true_indel_comp:
                            continue
                        if len(value[0]) > 1:
                            called_snp[key] = value
                        else:
                            snp_num += 1
                            if snp_num > expected_snp_num:
                                if aln_dis <= 0 and chr_dis < 600 and chr_dis > 200:
                                    called_snp[key] = value
                            elif snp_num > expected_snp_num - 5000:
                                if aln_dis <= 1 and chr_dis < 600 and chr_dis > 200:
                                    called_snp[key] = value
                            elif snp_num > expected_snp_num - 10000:
                                if aln_dis <= 1 and chr_dis < 600 and chr_dis > 200:
                                    called_snp[key] = value
                                elif aln_dis <= 2 and chr_dis < 575 and chr_dis > 225:
                                    called_snp[key] = value
                            elif snp_num > expected_snp_num - 15000:
                                if aln_dis <= 1 and chr_dis < 600 and chr_dis > 200:
                                    called_snp[key] = value
                                elif aln_dis <= 2 and chr_dis < 575 and chr_dis > 225:
                                    called_snp[key] = value
                                elif aln_dis <= 3 and chr_dis < 550 and chr_dis > 250:
                                    called_snp[key] = value
                            else:
                                called_snp[key] = value

                    print "# called snps after filter2", len(called_snp)

                    result_file.write("\t".join([prog_version, "%.0f" % (2.0*int(rn)*int(rl)/ref_len), str(confi)]) + "\t")
                    TP_KAKS, TP_NS = 0, 0
                    FP_KAKS, FP_NS = 0, 0
                    TP_KAKID, TP_NID = 0, 0
                    FP_KAKID, FP_NID = 0, 0
                    FP_S, FP_ID = 0, 0
                    for key, value in called_snp.iteritems():
                        if key in true_snp_comp or key in true_indel_comp:
                            if key in true_snp_comp:
                                if value[0] == true_snp_comp[key]:
                                    TP_KAKS += 1
                                else:
                                    FP_KAKS += 1
                            elif key in true_indel_comp:
                                if value[0] == true_indel_comp[key]:
                                    TP_KAKID += 1
                                else:
                                    FP_KAKID += 1
                        elif key in true_snp_none or key in true_indel_none:
                            if key in true_snp_none:
                                snp_num += 1
                                if value[0] == true_snp_none[key]:
                                    TP_NS += 1
                                else:
                                    FP_NS += 1
                            elif key in true_indel_none:
                                if value[0] == true_indel_none[key]:
                                    TP_NID += 1
                                else:
                                    FP_NID += 1
                        else:
                            if len(value[0]) == 1:
                                snp_num += 1
                                FP_S += 1
                            else:
                                FP_ID += 1

                    result_file.write(str(KAKS + NS) + "\t")
                    if TP_KAKS + FP_KAKS + TP_NS + FP_NS + FP_S != 0 and KAKS + NS != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKS + TP_NS)/float(TP_KAKS + TP_NS + FP_KAKS + FP_NS + FP_S)))
                        result_file.write("%.5f\t" % (float(TP_KAKS + TP_NS)/float(KAKS + NS)))
                    else:
                        result_file.write("\t\t")

                    result_file.write(str(FP_S) + "\t")
                    if TP_NS + FP_NS != 0 and NS != 0:
                        result_file.write("%.5f\t" % (float(TP_NS)/float(TP_NS + FP_NS)))
                        result_file.write("%.5f\t" % (float(TP_NS)/float(NS)))
                    else:
                        result_file.write("\t\t")
                    if TP_KAKS + FP_KAKS != 0 and KAKS != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKS)/float(TP_KAKS + FP_KAKS)))
                        result_file.write("%.5f\t" % (float(TP_KAKS)/float(KAKS)))
                    else:
                        result_file.write("\t\t")

                    result_file.write(str(KAKID + NID) + "\t")
                    if TP_KAKID + FP_KAKID + TP_NID + FP_NID + FP_ID != 0 and KAKID + NID != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKID + TP_NID)/float(TP_KAKID + TP_NID + FP_KAKID + FP_NID + FP_ID)))
                        result_file.write("%.5f\t" % (float(TP_KAKID + TP_NID)/float(KAKID + NID)))
                    else:
                        result_file.write("\t\t")

                    result_file.write(str(FP_ID) + "\t")
                    if TP_NID + FP_NID != 0 and NID != 0:
                        result_file.write("%.5f\t" % (float(TP_NID)/float(TP_NID + FP_NID)))
                        result_file.write("%.5f\t" % (float(TP_NID)/float(NID)))
                    else:
                        result_file.write("\t\t")
                    if TP_KAKID + FP_KAKID != 0 and KAKID != 0:
                        result_file.write("%.5f\t" % (float(TP_KAKID)/float(TP_KAKID + FP_KAKID)))
                        result_file.write("%.5f\t" % (float(TP_KAKID)/float(KAKID)))
                    else:
                        result_file.write("\t\t")

                    result_file.write(str(NS) + "\t" + str(TP_NS) + "\t" + str(FP_NS) + "\t")
                    result_file.write(str(KAKS) + "\t" + str(TP_KAKS) + "\t" + str(FP_KAKS) + "\t")

                    result_file.write(str(NID) + "\t" + str(TP_NID) + "\t" + str(FP_NID) + "\t")
                    result_file.write(str(KAKID) + "\t" + str(TP_KAKID) + "\t" + str(FP_KAKID) + "\t")

                    result_file.write(result_dn + "\t" + prefix_fn + "\t" + cpu_num + "\t" + str(ms) + "\t1\t")

                    mem_time_file = os.path.join(result_path, prefix_fn + ".snpcall." + str(cpu_num) + ".log")
                    with open(mem_time_file) as f:
                        for line in f:
                            tokens = line.strip().split("\t")
                            if "# of aligned reads" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                            if "# of no-aligned reads" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
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
                    with open(mem_time_file) as f:
                        for line in f:
                            tokens = line.strip().split("\t")
                            if "Input parameters" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                            if "Parameters" in tokens[0]:
                                result_file.write(tokens[1] + "\t")
                    result_file.write("\n")

    result_file.close()
