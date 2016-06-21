'''
Evaluate variant call results of Scalpel for real data
'''

import os
import sys
import bisect

if len(sys.argv) != 6:
    print "Usage: python eval_var_real_data_scalpel.py chr_name cov confi_S confi_I bed_file"
    exit(0)

chr_name = sys.argv[1]
cov = sys.argv[2]
confi_S = float(sys.argv[3])
confi_I = float(sys.argv[4])
bed_file = sys.argv[5]

itv1, itv2 = [], []
with open(bed_file) as f:
    for line in f:
        tokens = line.split("\t")
        itv1.append(int(tokens[1]))
        itv2.append(int(tokens[2]))
print "#intervals", len(itv1)

var_prof_file = os.path.join("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NISTIntegratedCalls.vcf")
var_prof = {}
TS, TI = 0, 0
with open(var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split("\t")
            
            #ignore het var
            if value[9].split(":")[0] != "1/1" and value[9].split(":")[0] != "1|1":
                continue
            
            #if value[0] == "chr" + chr_name:
            if value[0] == chr_name:
                #var_prof[int(value[1]) - 1] = value[3:5] + value[9:10]
                var_prof[int(value[1]) - 1] = value[3:5] + value[9].split(":")[0:1]
                if len(value[3]) == 1 and len(value[4]) == 1:
                    TS += 1
                else:
                    TI += 1
print "#var prof variants", TS, TI, TS + TI

dbsnp_var_prof = {}
dbsnp_var_prof_file = os.path.join("/backup2/nsvo/variant_calling/Human_data/refs/GRCh37_chr" + chr_name + "/TRIMMED.ALL.chr" + chr_name + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
#dbsnp_var_prof_file = "/data/nsvo/test_data/GRCh37_chr1/refs/ExAC.r0.3.1.sites.vep.vcf_chr1"
DS_KS, DS_KI = 0, 0
with open(dbsnp_var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split("\t")
            if value[0] == chr_name:
                dbsnp_var_prof[int(value[1]) - 1] = value[3:5]
                if len(value[3]) == 1 and len(value[4]) == 1:
                    DS_KS += 1
                else:
                    DS_KI += 1
print "#var dbsnp variants", DS_KS, DS_KI, DS_KS + DS_KI

result_path = os.path.join("/data/nsvo/test_data/GRCh37_chr" + chr_name + "/results/real_reads/NA12878")

var_call_file = os.path.join(result_path, "scalpel/" + cov + "/variants.indel.vcf")
var_call = {}
with open(var_call_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != '#':
            value = line.strip().split("\t")
            
            #ignore het var
            if value[9].split(":")[0] != "1/1" and value[9].split(":")[0] != "1|1":
                continue
            
            pos = int(value[1]) - 1
            if len(value[3]) == 1 and len(value[4]) == 1:
                #if float(value[5]) >= confi_S:
                    var_call[pos] = value[3:6] + value[9].split(":")[0:1]
            else:
                #if float(value[5]) >= confi_I:
                    var_call[pos] = value[3:6] + value[9].split(":")[0:1]
print "#called variants", len(var_call)

'''
cmp_tool_var, cmp_tool_var_all = {}, {}
confi_S_K, confi_I_K, confi_S_U, confi_I_U = 20, 20, 20, 20
cmp_tool_file = os.path.join(result_path, "ivc_0.8.1-3_exome", "ERR194147_shuf_50x_chr" + chr_name + ".ivc.vcf")
with open(cmp_tool_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split("\t")
            pos = int(value[1]) - 1
            i = bisect.bisect_right(itv1, pos)
            if i: # i.e., itv1[i-1] <= pos
                if pos > itv2[i-1]:
                    continue
            cmp_tool_var_all[pos] = line
            if pos in dbsnp_var_prof:
                if len(value[3]) == 1 and len(value[4]) == 1:
                    if float(value[5]) >= confi_S_U:
                        if float(value[12])/float(value[13]) >= 0.5 and int(value[12]) >= 5:
                            cmp_tool_var[pos] = value[3:6] + value[9:]
                else:
                    if float(value[5]) >= confi_S_U:
                        if float(value[12])/float(value[13]) >= 0.5 and int(value[12]) >= 5:
                            cmp_tool_var[pos] = value[3:6] + value[9:]
            else:
                if len(value[3]) == 1 and len(value[4]) == 1:
                    if float(value[5]) >= confi_S_U:
                        if float(value[12])/float(value[13]) >= 0.5 and int(value[12]) >= 5:
                            cmp_tool_var[pos] = value[3:6] + value[9:]
                else:
                    if float(value[5]) >= confi_I_U:
                        if float(value[12])/float(value[13]) >= 0.5 and int(value[12]) >= 5:
                            cmp_tool_var[pos] = value[3:6] + value[9:]
print "#cmp variants", len(cmp_tool_var), len(cmp_tool_var_all)
'''
'''
#cmp_tool_file = os.path.join(result_path, "samtools/ERR194147_shuf_" + cov + "_chr" + chr_name + ".bwa_sorted_RG_realign.norm.vcf")
cmp_tool_file = os.path.join(result_path, "gatk_hc/ERR194147_shuf_" + cov + "_chr" + chr_name + ".bwa_sorted_RG_realign.vcf")
cmp_tool_var, cmp_tool_var_all = {}, {}
with open(cmp_tool_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split("\t")
            cmp_tool_var_all[pos] = line
            
            if value[9].split(":")[0] != "1/1" and value[9].split(":")[0] != "1|1":
                continue
            
            if float(value[5]) >= 20:
                cmp_tool_var[int(value[1]) - 1] = value[3:6] + value[9].split(":")[0:1]
print "#cmp variants", len(cmp_tool_var)
'''
'''
analysis_dir = os.path.join(result_path, "scalpel/" + cov + "/cmp_to_ivc0.8.1_exome_" + str(confi_S_K) + "." + str(confi_I_K) + "." + str(confi_S_U) + "." + str(confi_I_U) + ".debug_" + cov + "_giab")
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

out_file11= open(os.path.join(analysis_dir, "ivc_miscall_tp_dbsnp_snp"), "w")
out_file21= open(os.path.join(analysis_dir, "ivc_notcall_tp_dbsnp_snp"), "w")
out_file31= open(os.path.join(analysis_dir, "ivc_miscall_tp_dbsnp_indel"), "w")
out_file41= open(os.path.join(analysis_dir, "ivc_notcall_tp_dbsnp_indel"), "w")

out_file22= open(os.path.join(analysis_dir, "ivc_lowqual_tp_snp"), "w")
out_file42= open(os.path.join(analysis_dir, "ivc_lowqual_tp_indel"), "w")

out_file91 = open(os.path.join(analysis_dir, "gatk_call_tp_na_dbsnp_snp"), "w")
out_file92 = open(os.path.join(analysis_dir, "gatk_call_tp_na_dbsnp_indel"), "w")
'''
KA_KS, NA_KS, NA_KS_DS, KA_KI, NA_KI, NA_KI_DI, NS, NS_DS, NI, NI_DI = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
KA_KS_DS, KA_KS_ST, KA_KI_DI, KA_KI_ST, NS_DS_ST, NS_ST, NI_DI_ST, NI_ST = 0, 0, 0, 0, 0, 0, 0, 0
for var_pos, var in var_call.iteritems():
    if var_pos in var_prof:
        if len(var[0]) == 1 and len(var[1]) == 1:
            if var[0] == var_prof[var_pos][0] and var[1] in set(var_prof[var_pos][1].split(",")):
                if var_pos in dbsnp_var_prof and var[0] == dbsnp_var_prof[var_pos][0] and var[1] in set(dbsnp_var_prof[var_pos][1].split(",")):
                    KA_KS_DS += 1
                else:
                    KA_KS += 1
                '''
                if var_pos in cmp_tool_var and set(var[0:2]) == set(cmp_tool_var[var_pos][0:2]):
                    KA_KS_ST += 1

                if var_pos in cmp_tool_var and set(var[0:2]) != set(cmp_tool_var[var_pos][0].split("|") + cmp_tool_var[var_pos][1].split("|")):
                    out_file1.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                    if var_pos in dbsnp_var_prof and var[0:2] == dbsnp_var_prof[var_pos]:
                        out_file11.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                if var_pos not in cmp_tool_var:
                    if var_pos not in cmp_tool_var_all:
                        out_file2.write(str(var_pos + 1) + "\t" + str(var) + "\n")
                    else:
                        out_file22.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var_all[var_pos]))
                    if var_pos in dbsnp_var_prof and var[0:2] == dbsnp_var_prof[var_pos]:
                        out_file21.write(str(var_pos + 1) + "\t" + str(var) + "\n")

                if var_pos in dbsnp_var_prof and var[0:2] != dbsnp_var_prof[var_pos]:
                    out_file91.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(dbsnp_var_prof[var_pos]) + "\n")
                '''
            else:
                if var_pos in dbsnp_var_prof and var[0] == dbsnp_var_prof[var_pos][0] and var[1] in set(dbsnp_var_prof[var_pos][1].split(",")):
                    NA_KS_DS += 1
                else:
                    NA_KS += 1
        else:
            if var[0] == var_prof[var_pos][0] and var[1] in set(var_prof[var_pos][1].split(",")):
                if var_pos in dbsnp_var_prof and var[0] == dbsnp_var_prof[var_pos][0] and var[1] in set(dbsnp_var_prof[var_pos][1].split(",")):
                    KA_KI_DI += 1
                else:
                    KA_KI += 1
                '''
                if var_pos in cmp_tool_var and set(var[0:2]) == set(cmp_tool_var[var_pos][0:2]):
                    KA_KI_ST += 1

                if var_pos in cmp_tool_var and set(var[0:2]) != set(cmp_tool_var[var_pos][0].split("|") + cmp_tool_var[var_pos][1].split("|")):
                    out_file3.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                    if var_pos in dbsnp_var_prof and var[0:2] == dbsnp_var_prof[var_pos]:
                        out_file31.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                if var_pos not in cmp_tool_var:
                    if var_pos not in cmp_tool_var_all:
                        out_file4.write(str(var_pos + 1) + "\t" + str(var) + "\n")
                    else:
                        out_file42.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var_all[var_pos]))
                    if var_pos in dbsnp_var_prof and var[0:2] == dbsnp_var_prof[var_pos]:
                        out_file41.write(str(var_pos + 1) + "\t" + str(var) + "\n")

                if var_pos in dbsnp_var_prof and var[0:2] != dbsnp_var_prof[var_pos]:
                    out_file92.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(dbsnp_var_prof[var_pos]) + "\n")
                '''
            else:
                if var_pos in dbsnp_var_prof and var[0] == dbsnp_var_prof[var_pos][0] and var[1] in set(dbsnp_var_prof[var_pos][1].split(",")):
                    NA_KI_DI += 1
                else:
                    NA_KI += 1
    else:
        if len(var[0]) == 1 and len(var[1]) == 1:
            if var_pos in dbsnp_var_prof and var[0] == dbsnp_var_prof[var_pos][0] and var[1] in set(dbsnp_var_prof[var_pos][1].split(",")):
                NS_DS += 1
                '''
                if var_pos in cmp_tool_var and set(var[0:2]) == set(cmp_tool_var[var_pos][0:2]):
                    NS_DS_ST += 1
                '''
            else:
                NS += 1
                '''
                if var_pos in cmp_tool_var and set(var[0:2]) == set(cmp_tool_var[var_pos][0:2]):
                    NS_ST += 1
                if var_pos in cmp_tool_var:
                    out_file5.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                else:
                    out_file6.write(str(var_pos + 1) + "\t" + str(var) + "\n")
                '''
        else:
            if var_pos in dbsnp_var_prof and var[0] == dbsnp_var_prof[var_pos][0] and var[1] in set(dbsnp_var_prof[var_pos][1].split(",")):
                NI_DI += 1
                '''
                if var_pos in cmp_tool_var and set(var[0:2]) == set(cmp_tool_var[var_pos][0:2]):
                    NI_DI_ST += 1
                '''
            else:
                NI += 1
                '''
                if var_pos in cmp_tool_var and set(var[0:2]) == set(cmp_tool_var[var_pos][0:2]):
                    NI_ST += 1
                if var_pos in cmp_tool_var:
                    out_file7.write(str(var_pos + 1) + "\t" + str(var) + "\t" + str(cmp_tool_var[var_pos]) + "\n")
                else:
                    out_file8.write(str(var_pos + 1) + "\t" + str(var) + "\n")
                '''
'''
out_file1.close()
out_file2.close()
out_file3.close()
out_file4.close()

out_file5.close()
out_file6.close()
out_file7.close()
out_file8.close()

out_file22.close()
out_file42.close()

out_file91.close()
out_file92.close()
'''
CS = KA_KS + KA_KS_DS + NA_KS + NA_KS_DS + NS_DS + NS
CI = KA_KI + KA_KI_DI + NA_KI + NA_KI_DI + NI_DI + NI

result_file_path = os.path.join(result_path, "scalpel/" + cov + "/variants.indel.vcf.homvar.prec_rec.GIAB." + str(confi_S) + "." + str(confi_I) + ".txt")
result_file = open(result_file_path, "w")

header = ["Alg", "Cov", "Qual", "TP_U_S", "TP_K_S", "TL_FP_U_S", "TL_FP_K_S", "TP_U_I", "TP_K_I", "TL_FP_U_I", "TL_FP_K_I", "FL_FP_U_S", "FL_FP_K_S", "FL_FP_U_I", "FL_FP_K_I", "TS", "TI", "TS+TI", "CS", "CI", "CS+CI", "PS", "RS", "PI", "RI"]

result_file.write("\t".join(header) + "\n")
result_file.write("Scalpel\t" + cov + "\t" + str(confi_S) + "," + str(confi_I) + "\t")
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t" % (KA_KS, KA_KS_DS, NA_KS, NA_KS_DS, KA_KI, KA_KI_DI, NA_KI, NA_KI_DI, NS, NS_DS, NI, NI_DI))
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5f\t%.5f\t%.5f\t%.5f" % (TS, TI, TS + TI, CS, CI, CS+CI, KA_KS/float(CS+1), KA_KS/float(TS), KA_KI/float(CI+1), KA_KI/float(TI)))
result_file.close()

'''
header = ["Alg", "Cov", "Qual", "TP_U_S", "TP_K_S", "TL_FP_U_S", "TL_FP_K_S", "TP_U_I", "TP_K_I", "TL_FP_U_I", "TL_FP_K_I", "FL_FP_U_S", "FL_FP_K_S", "FL_FP_U_I", "FL_FP_K_I", "TS", "TI", "TS+TI", "CS", "CI", "CS+CI", "PS", "RS", "PI", "RI"]
result_file.write("\t".join(header) + "\n")
result_file.write("Scalpel\t" + cov + "\t" + str(confi_S) + "," + str(confi_I) + "\t")
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t" % (KA_KS, KA_KS_DS, KA_KS_ST, NA_KS, NA_KS_DS, KA_KI, KA_KI_DI, KA_KI_ST, NA_KI, NA_KI_DI, NS_DS, NS_DS_ST, NI_DI, NI_DI_ST, NS, NS_ST, NI, NI_ST))
result_file.write("%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5d\t%.5f\t%.5f\t%.5f\t%.5f" % (TS, TI, TS + TI, CS, CI, CS+CI, KA_KS/float(S), KA_KS/float(TS), KA_KI/float(CI), KA_KI/float(TI)))
result_file.close()
'''
print "Check results at:", result_file_path
