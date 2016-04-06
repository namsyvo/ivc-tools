'''
Evaluate variant call results for real data
'''

import os
import sys
import json
from datetime import datetime

if len(sys.argv) != 4:
    print "Usage: python eval_var_af_sid_mutant_diff_ref.py confi_S confi_I chr_name"
    exit(0)

confi_S = float(sys.argv[1])
confi_I = float(sys.argv[2])
chr_name = sys.argv[3]

#var_prof_file = os.path.join("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NA12878.vcf")
var_prof_file = os.path.join("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NISTIntegratedCalls.vcf")
var_prof = {}
KS, KI = 0, 0
with open(var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            
            #ignore het var
            if value[9].split(":")[0] != "1/1" and value[9].split(":")[0] != "1|1":
                continue
            
            #if value[0] == "chr" + chr_name:
            if value[0] == chr_name:
                #var_prof[int(value[1]) - 1] = value[3:5] + value[9:10]
                var_prof[int(value[1]) - 1] = value[3:5] + value[9].split(":")[0:1]
                if len(value[3]) == 1 and len(value[4]) == 1:
                    KS += 1
                else:
                    KI += 1
print "#var prof variants", KS, KI, KS + KI

dbsnp_var_prof = {}
dbsnp_var_prof_file = os.path.join("/backup2/nsvo/variant_calling/Human_data/refs/GRCh37_chr" + chr_name + "/TRIMMED.ALL.chr" + chr_name + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
DS_KS, DS_KI = 0, 0
with open(dbsnp_var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            if value[0] == chr_name:
                dbsnp_var_prof[int(value[1]) - 1] = value[3:5]
                if len(value[3]) == 1 and len(value[4]) == 1:
                    DS_KS += 1
                else:
                    DS_KI += 1
print "#var dbsnp variants", DS_KS, DS_KI, DS_KS + DS_KI

result_path = os.path.join("/data/nsvo/test_data/GRCh37_chr" + chr_name + "/results/real_reads/NA12878")

var_call_file = os.path.join(result_path, "gatk_hc/ERR194147_sorted_chr" + chr_name + ".bwa_sorted_RG_realign.vcf")
#var_call_file = os.path.join(result_path, "samtools/ERR194147_sorted_chr" + chr_name + ".bwa_sorted_RG_realign.norm.vcf")
var_call = {}
with open(var_call_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split()
            
            #ignore het var
            if value[9].split(":")[0] != "1/1" and value[9].split(":")[0] != "1|1":
                continue
            
            pos = int(value[1]) - 1
            if len(value[3]) == 1 and len(value[4]) == 1:
                if float(value[5]) >= confi_S:
                    var_call[pos] = value[3:5] + value[9].split(":")[0:1]
            else:
                if float(value[5]) >= confi_I:
                    var_call[pos] = value[3:5] + value[9].split(":")[0:1]
print "#called variants", len(var_call)

cmp_tool_file = os.path.join(result_path, "ivc_0.9.0/ERR194147_sorted_chr" + chr_name + ".ivc_0.9.0.vcf")
cmp_tool_var = {}
with open(cmp_tool_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            if value[5] == "NaN" or float(value[5]) >= 20:
                cmp_tool_var[int(value[1]) - 1] = value[3:5]
print "#cmp variants", len(cmp_tool_var)

analysis_dir = os.path.join(result_path, "gatk_hc/cmp_to_ivc_nist")
#analysis_dir = os.path.join(result_path, "samtools/cmp_to_ivc_nist")
if not os.path.exists(analysis_dir):
    os.mkdir(analysis_dir)

out_file1 = open(os.path.join(analysis_dir, "ivc_miscall_tp_snp"), "w")
out_file2 = open(os.path.join(analysis_dir, "ivc_notcall_tp_snp"), "w")
out_file3 = open(os.path.join(analysis_dir, "ivc_miscall_tp_indel"), "w")
out_file4 = open(os.path.join(analysis_dir, "ivc_notcall_tp_indel"), "w")

out_file5 = open(os.path.join(analysis_dir, "ivc_call_ns"), "w")
out_file6 = open(os.path.join(analysis_dir, "ivc_notcall_ns"), "w")
out_file7 = open(os.path.join(analysis_dir, "ivc_call_ni"), "w")
out_file8 = open(os.path.join(analysis_dir, "ivc_notcall_ni"), "w")

KA_KS, NA_KS, KA_KI, NA_KI, NS, NI = 0, 0, 0, 0, 0, 0
NS_DS, NI_DI = 0, 0
for var_pos, var in var_call.iteritems():
    if var_pos in var_prof:
        if len(var[0]) == 1 and len(var[1]) == 1:
            if var[0:2] == var_prof[var_pos][0:2]:
                KA_KS += 1
                if var_pos in cmp_tool_var and set(var_prof[var_pos][0:2]) != set(cmp_tool_var[var_pos][0].split("|") + cmp_tool_var[var_pos][1].split("|")):
                    out_file1.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                elif var_pos not in cmp_tool_var:
                    out_file2.write(str(var_pos + 1) + "\t" + str(var) + "\n")
            else:
                NA_KS += 1
        else:
            if var[0:2] == var_prof[var_pos][0:2]:
                KA_KI += 1
                if var_pos in cmp_tool_var and set(var_prof[var_pos][0:2]) != set(cmp_tool_var[var_pos][0].split("|") + cmp_tool_var[var_pos][1].split("|")):
                    out_file3.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                elif var_pos not in cmp_tool_var:
                    out_file4.write(str(var_pos + 1) + "\t" + str(var) + "\n")
            else:
                NA_KI += 1
    else:
        if len(var[0]) == 1 and len(var[1]) == 1:
            NS += 1
            if var_pos in dbsnp_var_prof:
                if var[0:2] == dbsnp_var_prof[var_pos]:
                    NS_DS += 1
            else:
                if var_pos in cmp_tool_var:
                    out_file5.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                else:
                    out_file6.write(str(var_pos + 1) + "\t" + str(var) + "\n")
        else:
            NI += 1
            if var_pos in dbsnp_var_prof:
                if var[0:2] == dbsnp_var_prof[var_pos]:
                    NI_DI += 1
            else:
                if var_pos in cmp_tool_var:
                    out_file7.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                else:
                    out_file8.write(str(var_pos + 1) + "\t" + str(var) + "\n")

out_file1.close()
out_file2.close()
out_file3.close()
out_file4.close()

out_file5.close()
out_file6.close()
out_file7.close()
out_file8.close()

S = KA_KS + NA_KS + NS
I = KA_KI + NA_KI + NI

#result_file_path = os.path.join(result_path, "gatk_hc/ERR194147_sorted_chr" + chr_name + ".bwa_sorted_RG_realign.vcf.prec_rec.Platinum." + str(confi_S) + "." + str(confi_I) + ".txt")
result_file_path = os.path.join(result_path, "gatk_hc/ERR194147_sorted_chr" + chr_name + ".bwa_sorted_RG_realign.vcf.prec_rec.GIAB." + str(confi_S) + "." + str(confi_I) + ".txt")
#result_file_path = os.path.join(result_path, "samtools/ERR194147_sorted_chr" + chr_name + ".bwa_sorted_RG_realign.norm.vcf.prec_rec.Platinum." + str(confi_S) + "." + str(confi_I) + ".txt")
#result_file_path = os.path.join(result_path, "samtools/ERR194147_sorted_chr" + chr_name + ".bwa_sorted_RG_realign.norm.vcf.prec_rec.GIAB." + str(confi_S) + "." + str(confi_I) + ".txt")
result_file = open(result_file_path, "w")

header = ["Alg", "Cov", "Qual", "KA_KS", "NA_KS", "KA_KI", "NA_KI", "NS", "NI", "NS_DS", "NI_DI", "KS", "KI", "KS+KI", "CS", "CI", "CS+CI", "PS", "RS", "PI", "RI"]

result_file.write("\t".join(header) + "\n")
result_file.write("GATK-HC\t50\t" + str(confi_S) + "," + str(confi_I) + "\t")
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t" % (KA_KS, NA_KS, KA_KI, NA_KI, NS, NI, NS_DS, NI_DI, KS, KI, KS + KI, S, I, S+I))
result_file.write("%.5f\t%.5f\t%.5f\t%.5f" % (KA_KS/float(S), KA_KS/float(KS), KA_KI/float(I), KA_KI/float(KI)))
result_file.close()
print "Check results at:", result_file_path
