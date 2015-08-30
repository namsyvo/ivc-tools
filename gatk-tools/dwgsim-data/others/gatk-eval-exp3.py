'''
Evaluate GATK performance
Consider two categories of SNP/Indels based on reference SNP profile
Do not evaluate SNPs/Indels which are same as Ref
Consider only SNPs and Insertions (no Deletions)
Usage: python gatk-eval-exp3.py config-dwgsim-chr1.txt 20
'''

import os
import sys
from datetime import datetime

config_file = sys.argv[1]
f=open(config_file)
prog_path = f.readline().strip()
data_path = f.readline().strip()
genome_fn = f.readline().strip()
snp_fn = f.readline().strip()
read_fn = f.readline().strip()
f.close()

confi = float(sys.argv[2])

ref_len = 249250621
read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]]

known_true_snp, known_true_indel, unknown_true_snp, unknown_true_indel, var_prof = {}, {}, {}, {}, {}

ref_path = os.path.join(data_path,  "refs", "af_mutant")
variant_comp_file = os.path.join(ref_path, "known_var_0.7.txt")
variant_none_file = os.path.join(ref_path, "unknown_var_0.7.txt")
variant_prof_file = os.path.join(data_path, "refs", "TRIMMED.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.diffcontigname.vcf")

with open(variant_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            var_prof[int(value[1]) - 1] = value[3:5]

KAKS, NS = 0, 0
KAKID, NID = 0, 0

with open(variant_comp_file) as f:
    for line in f.readlines():
        if line.strip():
            value = line.strip().split()
            pos, val = int(value[0]), value[1]
            known_var = var_prof[pos]
            if len(known_var[0]) == 1 and len(known_var[1]) == 1:
                known_true_snp[pos] = val
                if val != known_var[0]:
                    KAKS += 1
            else:
                known_true_indel[pos] = val
                if val != known_var[0]:
                    KAKID += 1

with open(variant_none_file) as f:
    for line in f.readlines():
        if line.strip():
            value = line.strip().split()
            pos, val = int(value[0]), value[1]
            unknown_var = var_prof[pos]
            if len(unknown_var[0]) == 1 and len(unknown_var[1]) == 1:
                unknown_true_snp[pos] = val
                if val != unknown_var[0]:
                    NS += 1
            else:
                unknown_true_indel[pos] = val
                if val != unknown_var[0]:
                    NID += 1

result_path = os.path.join(data_path, "results/sim-reads/af_mutant_dwgsim/gatk")
result_file_path = result_path + "/" + read_fn + "_" + str(read_lens[0]) + "." + str(seq_errs[0]) + ".prec-rec-time-mem." + str(confi) + ".txt"
result_file = open(result_file_path, "w")

header = ["Alg", "cov", "qual", "TP-S", "FP-S", "FP-S-N", "TP-S-U", "FP-S-U", "TP-S-K", "FP-S-K", \
            "TP-I", "FP-I", "FP-I-N", "TP-I-U", "FP-I-U", "TP-I-K", "FP-I-K", \
            "P-S", "R-S", "P-S-U", "R-S-U", "P-S-K", "R-S-K", "P-I", "R-I", "P-I-U", "R-I-U", "P-I-K", "R-I-K", \
            "S", "S-U", "S-K", "CS", "CS-S", "I", "I-U", "I-K", "CI", "CI-I", "run", "read"]

result_file.write("\t".join(header))
result_file.write("\n")

time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "-")

for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            prefix_fn = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
            called_snp_file = result_path + "/" + prefix_fn + ".bwa.vcf"
            snp = {}
            f = open(called_snp_file)
            for line in f.readlines():
                if line[0] != '#':
                    value = line.strip().split()
                    if float(value[5]) >= confi:
                        snp[int(value[1]) - 1] = value[4]

            TP_KAKS, TP_NS = 0, 0
            FP_KAKS, FP_NS = 0, 0
            TP_KAKID, TP_NID = 0, 0
            FP_KAKID, FP_NID = 0, 0
            FP_S, FP_ID = 0, 0
            for key, value in snp.iteritems():
                if key in known_true_snp or key in known_true_indel:
                    if key in known_true_snp:
                        if value == known_true_snp[key]:
                            TP_KAKS += 1
                        else:
                            FP_KAKS += 1
                    elif key in known_true_indel:
                        if value == known_true_indel[key]:
                            TP_KAKID += 1
                        else:
                            FP_KAKID += 1
                elif key in unknown_true_snp or key in unknown_true_indel:
                    if key in unknown_true_snp:
                        if value == unknown_true_snp[key]:
                            TP_NS += 1
                        else:
                            FP_NS += 1
                    elif key in unknown_true_indel:
                        if value == unknown_true_indel[key]:
                            TP_NID += 1
                        else:
                            FP_NID += 1
                else:
                        if len(value) == 1:
                            FP_S += 1
                        else:
                            FP_ID += 1

            '''
            nums_header = ["Alg", "cov", "qual", "TP-S", "FP-S", "FP-S-N", "TP-S-U", "FP-S-U", "TP-S-K", "FP-S-K", \
                        "TP-I", "FP-I", "FP-I-N", "TP-I-U", "FP-I-U", "TP-I-K", "FP-I-K"]
            '''
            result_file.write("\t".join(["GATK-3.1-1", "%.0f" % (2.0*int(rn)*int(rl)/ref_len), str(confi)]) + "\t")

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
            para_header = ["S", "S-U", "S-K", "CS", "CS\S", "I", "I-U", "I-K", "CI", "CI\I", "run", "read"]
            '''
            #S, "S-U", "S-K", CS, CS-S
            result_file.write(str(KAKS + NS) + "\t" + str(NS) + "\t" + str(KAKS) + "\t")
            result_file.write("%.5d\t" % (FP_KAKS + FP_NS + FP_S + TP_KAKS + TP_NS))
            result_file.write("%.5d\t" % ((FP_KAKS + FP_NS + FP_S + TP_KAKS + TP_NS) - (KAKS + NS)))

            #"I", "I-U", "I-K", "CI", "CI-I"
            result_file.write(str(KAKID + NID) + "\t" + str(NID) + "\t" + str(KAKID) + "\t")
            result_file.write("%.5d\t" % (FP_KAKID + FP_NID + FP_ID + TP_KAKID + TP_NID))
            result_file.write("%.5d\t" % ((FP_KAKID + FP_NID + FP_ID + TP_KAKID + TP_NID) - (KAKID + NID)))

            #"run", "read"
            result_file.write("GATK-" + time_stamp + "\t" + prefix_fn + "\t")
            result_file.write("\n")

result_file.close()
