'''
Evaluation script for variant calling
Usage: python gatk-eval-exp.py  <snp call file>  <quality threshold>  <true snp file>  <result file>
E.g.: python gatk-eval-exp.py snp_call.vcf 20 true_snp.vcf eval.txt
'''

import os
import sys

snp = {}
f = open(sys.argv[1])
confi = float(sys.argv[2])
for line in f.readlines():
    if line[0] != '#':
        value = line.strip().split()
        if float(value[5]) >= confi:
            snp[int(value[1])] = value[4]

true_snp, true_indel = {}, {}
f = open(sys.argv[3])
for line in f.readlines():
    if line[0] != '#':
        value = line.strip().split()
        if len(value[4]) == 1 and value[4] != ".":
            true_snp[int(value[1])] = value[4]
        else:
            true_indel[int(value[1])] = value[4]

result_file = open(sys.argv[4], "w")

header = ["Confi", "Prec", "Rec", "Prec-SNP", "Rec-SNP", "Prec-Indel", "Rec-Indel", 
          "#SNP", "#Indel", "TP-SNP", "FP-SNP", "TP-Indel", "FP-Indel", "FP_NS", "FP_NID"]

result_file.write("\t".join(header) + "\n")
result_file.write(str(confi) + "\t")

KS, KID = len(true_snp), len(true_indel)
TP_KS, FP_KS, TP_KID, FP_KID, FP_NS, FP_NID = 0, 0, 0, 0, 0, 0

for key, value in snp.iteritems():
    if key in true_snp or key in true_indel:
        if key in true_snp:
            if value == true_snp[key]:
                TP_KS += 1
            else:
                FP_KS += 1
        elif key in true_indel:
            if value == true_indel[key]:
                TP_KID += 1
            else:
                FP_KID += 1
    else:
        if len(value) == 1:
            FP_NS += 1
        else:
            FP_NID += 1

result_file.write("%.5f\t" % (float(TP_KS + TP_KID)/float(TP_KS + FP_KS + FP_NS + TP_KID + FP_KID + FP_NID)))
result_file.write("%.5f\t" % (float(TP_KS + TP_KID)/float(KS + KID)))

if TP_KS + FP_KS + FP_NS != 0:
    result_file.write("%.5f\t" % (float(TP_KS)/float(TP_KS + FP_KS + FP_NS)))
    result_file.write("%.5f\t" % (float(TP_KS)/float(KS)))

if TP_KID + FP_KID + FP_NID != 0:
    result_file.write("%.5f\t" % (float(TP_KID)/float(TP_KID + FP_KID + FP_NID)))
    result_file.write("%.5f\t" % (float(TP_KID)/float(KID)))

result_file.write(str(KS) + "\t" + str(KID) + "\t")
result_file.write(str(TP_KS) + "\t" + str(FP_KS) + "\t")
result_file.write(str(TP_KID) + "\t" + str(FP_KID) + "\t")
result_file.write(str(FP_NS) + "\t" + str(FP_NID) + "\n")

result_file.close()
