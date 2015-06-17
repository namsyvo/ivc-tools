'''
Evaluate GATK performance
Consider two categories of SNP/Indels based on mutant-specific SNP profile
Evaluate SNPs/Indels which are both same as or different from Ref
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

ref_path = data_path + "/refs/mutant_af"
genome_file = ref_path + "/" + genome_fn

true_snp_comp, true_indel_comp, true_snp_part, true_indel_part, true_snp_none, true_indel_none = {}, {}, {}, {}, {}, {}

variant_comp_file = os.path.join(ref_path, "variant_comp_0.70.txt")
variant_none_file = os.path.join(ref_path, "variant_none_0.70.txt")

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

result_path = os.path.join(data_path, "results/sim-reads/mutant_af_dwgsim/gatk")
result_file_path = result_path + "/" + read_fn + "_" + str(read_lens[0]) + "." + str(seq_errs[0]) + ".prec-rec-time-mem." + str(confi) + ".txt"
result_file = open(result_file_path, "w")

header = ["Alg", "cov", "qual", \
            "S", "P-S", "R-S", "FP-S@Other", "P-S@None", "R-S@None", "P-S@Comp", "R-S@Comp", \
            "I", "P-I", "R-I", "FP-I@Other", "P-I@None", "R-I@None", "P-I@Comp", "R-I@Comp", \
            "S@None", "TP-S@None", "FP-S@None", "S@Comp", "TP-S@Comp", "FP-S@Comp", \
            "I@None", "TP-I@None", "FP-I@None", "I@Comp", "TP-I@Comp", "FP-I@Comp", \
            "run", "read", "proc", "timeI", "memI", "timeC", "memC"]

result_file.write("\t".join(header))
result_file.write("\n")

time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "-")

for rl in read_lens:
    for err in seq_errs:
        for rn in read_nums:
            result_file.write("GATK\t" + "%.0f\t" % (2.0*int(rn)*int(rl)/ref_len) + str(confi) + "\t")

            called_snp_file = result_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.vcf"
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

            result_file.write(time_stamp + "\t" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + "\t1\t0\t0\t0\t0\n")

result_file.close()
