'''
Extract TP, FP, and FN from SNP call results
Usage: python get_tpfpfn_var_diff_ref_af_sid_mutant.py config_file confi coverage_num result_dir_name
    E.g.: python get_tpfpfn_var_diff_ref_af_sid_mutant.py config-chr1-test.json 1.14 5 called_var_dir
'''
import os
import sys
import json

if len(sys.argv) != 4:
    print "Usage: python get_tpfpfn.py config_file confi_num coverage_num"
    exit(0)

config_file = sys.argv[1]
f=open(config_file)
prog_path = f.readline().strip()
data_path = f.readline().strip()
genome_fn = f.readline().strip()
dbsnp_fn = f.readline().strip()
read_fn = f.readline().strip()
f.close()

confi = float(sys.argv[2])
cov_num = sys.argv[3]

ref_len = 249250621
ref_para = ['0.70']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
max_snum = [2**i for i in range(3, 14)]
read_nums = []
if cov_num == "all":
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25]]
else:
    read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [int(cov_num)]]

var_prof = {}
var_prof_file = os.path.join(data_path, "refs/TRIMMED.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.diffcontigname.vcf")
with open(var_prof_file) as f:
    for line in f.readlines():
        if line.strip() and line[0] != "#":
            value = line.strip().split()
            var_prof[int(value[1]) - 1] = value[3:5]

for para in ref_para:
    true_known_snp, true_known_indel, true_unknown_snp, true_unknown_indel,  = {}, {}, {}, {}
    known_var_file = os.path.join(data_path, "refs/af_sid_mutant/known_var_" + para + ".txt")
    unknown_var_file = os.path.join(data_path, "refs/af_sid_mutant/unknown_var_" + para + ".txt")

    with open(known_var_file) as f:
        for line in f.readlines():
            if line.strip() and line[0] != '#':
                value = line.strip().split()
                pos, known_var = int(value[0]), value[1:]
                if len(var_prof[pos][0]) == 1 and len(var_prof[pos][1]) == 1:
                    true_known_snp[pos] = known_var
                else:
                    true_known_indel[pos] = known_var

    with open(unknown_var_file) as f:
        for line in f.readlines():
            if line.strip() and line[0] != '#':
                value = line.strip().split()
                pos, unknown_var = int(value[0]), value[1:]
                if len(var_prof[pos][0]) == 1 and len(var_prof[pos][1]) == 1:
                    true_unknown_snp[pos] = unknown_var
                else:
                    true_unknown_indel[pos] = unknown_var

    result_dn = os.path.join(data_path, "results/sim-reads/af_sid_mutant_dwgsim/gatk")

    fpfntp_info_path = os.path.join(result_dn, "fpfntp_info")
    if not os.path.exists(fpfntp_info_path):
        os.makedirs(fpfntp_info_path)

    header = "pos\tref\talt\tqual\tfilter\tinfo\tformat\tDepristo\n"
    fp_header = "pos\ttrue_ref\ttrue_alt\tref\talt\tqual\tfilter\tinfo\tformat\tDepristo\n"
    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                fn_part = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
                file_prefix = fpfntp_info_path + "/" + fn_part

                tp_snp_known_file = open(file_prefix + ".tp_snp_known." + str(confi) + ".txt", "w")
                tp_snp_known_file.write(header)
                fp_snp_known_file = open(file_prefix + ".fp_snp_known." + str(confi) + ".txt", "w")
                fp_snp_known_file.write(fp_header)
                tp_indel_known_file = open(file_prefix + ".tp_indel_known." + str(confi) + ".txt", "w")
                tp_indel_known_file.write(header)
                fp_indel_known_file = open(file_prefix + ".fp_indel_known." + str(confi) + ".txt", "w")
                fp_indel_known_file.write(fp_header)

                tp_snp_unknown_file = open(file_prefix + ".tp_snp_unknown." + str(confi) + ".txt", "w")
                tp_snp_unknown_file.write(header)
                fp_snp_unknown_file = open(file_prefix + ".fp_snp_unknown." + str(confi) + ".txt", "w")
                fp_snp_unknown_file.write(fp_header)

                tp_indel_unknown_file = open(file_prefix + ".tp_indel_unknown." + str(confi) + ".txt", "w")
                tp_indel_unknown_file.write(header)
                fp_indel_unknown_file = open(file_prefix + ".fp_indel_unknown." + str(confi) + ".txt", "w")
                fp_indel_unknown_file.write(fp_header)

                fp_snp_none_file = open(file_prefix + ".fp_snp_none." + str(confi) + ".txt", "w")
                fp_snp_none_file.write(header)
                fp_indel_none_file = open(file_prefix + ".fp_indel_none." + str(confi) + ".txt", "w")
                fp_indel_none_file.write(header)

                fn_snp_known_file = open(file_prefix + ".fn_snp_known." + str(confi) + ".txt", "w")
                fn_snp_known_file.write("pos\tsnp\n")
                fn_snp_unknown_file = open(file_prefix + ".fn_snp_unknown." + str(confi) + ".txt", "w")
                fn_snp_unknown_file.write("pos\tsnp\n")

                fn_indel_known_file = open(file_prefix + ".fn_indel_known." + str(confi) + ".txt", "w")
                fn_indel_known_file.write("pos\tsnp\n")
                fn_indel_unknown_file = open(file_prefix + ".fn_indel_unknown." + str(confi) + ".txt", "w")
                fn_indel_unknown_file.write("pos\tsnp\n")

                fn_snp_known_lowqual_file = open(file_prefix + ".fn_snp_known_lowqual." + str(confi) + ".txt", "w")
                fn_snp_unknown_lowqual_file = open(file_prefix + ".fn_snp_unknown_lowqual." + str(confi) + ".txt", "w")

                fn_indel_known_lowqual_file = open(file_prefix + ".fn_indel_known_lowqual." + str(confi) + ".txt", "w")
                fn_indel_unknown_lowqual_file = open(file_prefix + ".fn_indel_unknown_lowqual." + str(confi) + ".txt", "w")

                fp_low_qual_file = open(file_prefix + ".fp_low_qual." + str(confi) + ".txt", "w")

                called_var, low_qual_snp = {}, {}
                called_var_file = os.path.join(result_dn, fn_part + ".bwa.vcf")
                f = open(called_var_file)
                for line in f.readlines():
                    value = line.strip().split()
                    if value[0][0] != '#' and float(value[5]) >= confi:
                        pos = int(value[1]) - 1
                        var = value[3:5]
                        called_var[pos] = var
                        var_call_info = "\t".join(value[3:]) + "\n"
                        if pos in true_known_snp or pos in true_known_indel:
                            if value[3] != value[4]:
                                if pos in true_known_snp:
                                    if var == true_known_snp[pos]:
                                        tp_snp_known_file.write(str(pos) + "\t" + var_call_info)
                                    else:
                                        fp_snp_known_file.write(str(pos) + "\t" + true_known_snp[pos][0] + "\t" + true_known_snp[pos][1] + "\t" + var_call_info)
                                elif pos in true_known_indel:
                                    if var == true_known_indel[pos]:
                                        tp_indel_known_file.write(str(pos) + "\t" + var_call_info)
                                    else:
                                        fp_indel_known_file.write(str(pos) + "\t" + true_known_indel[pos][0] + "\t" + true_known_indel[pos][1] + "\t" + var_call_info)
                        elif pos in true_unknown_snp or pos in true_unknown_indel:
                            if pos in true_unknown_snp:
                                if var == true_unknown_snp[pos]:
                                    tp_snp_unknown_file.write(str(pos) + "\t" + var_call_info)
                                else:
                                    fp_snp_unknown_file.write(str(pos) + "\t" + true_unknown_snp[pos][0] + "\t" + true_unknown_snp[pos][1] + "\t" + var_call_info)
                            elif pos in true_unknown_indel:
                                if var == true_unknown_indel[pos]:
                                    tp_indel_unknown_file.write(str(pos) + "\t" + var_call_info)
                                else:
                                    fp_indel_unknown_file.write(str(pos) + "\t" + true_unknown_indel[pos][0] + "\t" + true_unknown_indel[pos][1] + "\t" + var_call_info)
                        else:
                            if len(value[3]) == 1 and len(value[4]) == 1:
                                fp_snp_none_file.write(str(pos) + "\t" + var_call_info)
                            else:
                                fp_indel_none_file.write(str(pos) + "\t" + var_call_info)
                    else:
                        low_qual_snp[pos] = value

                f.close()
                print "GATK", len(called_var)

                #Print FN info
                for pos, value in true_known_snp.iteritems():
                    if pos not in called_var and value != var_prof[pos]:
                        fn_snp_known_file.write(str(pos) + "\t" + "\t".join(value) + "\n")
                        if pos in low_qual_snp:
                            fn_snp_known_lowqual_file.write(str(pos) + "\t" + "\t".join(low_qual_snp[pos]) + "\n")

                for pos, value in true_unknown_snp.iteritems():
                    if pos not in called_var:
                        fn_snp_unknown_file.write(str(pos) + "\t" + "\t".join(value) + "\n")
                        if pos in low_qual_snp:
                            fn_snp_unknown_lowqual_file.write(str(pos) + "\t" + "\t".join(low_qual_snp[pos]) + "\n")

                for pos, value in true_known_indel.iteritems():
                    if pos not in called_var and value != var_prof[pos]:
                        fn_indel_known_file.write(str(pos) + "\t" + "\t".join(value) + "\n")
                        if pos in low_qual_snp:
                            fn_indel_known_lowqual_file.write(str(pos) + "\t" + "\t".join(low_qual_snp[pos]) + "\n")

                for pos, value in true_unknown_indel.iteritems():
                    if pos not in called_var:
                        fn_indel_unknown_file.write(str(pos) + "\t" + "\t".join(value) + "\n")
                        if pos in low_qual_snp:
                            fn_indel_unknown_lowqual_file.write(str(pos) + "\t" + "\t".join(low_qual_snp[pos]) + "\n")

                for pos, value in low_qual_snp.iteritems():
                    if pos not in true_known_snp and pos not in true_unknown_snp \
                        and pos not in true_known_indel and pos not in true_unknown_indel and pos in var_prof and value != var_prof[pos]:
                        fp_low_qual_file.write(str(pos) + "\t" + "\t".join(value) + "\n")

                tp_snp_known_file.close()
                fp_snp_known_file.close()
                tp_snp_unknown_file.close()
                fp_snp_unknown_file.close()
                
                tp_indel_known_file.close()
                fp_indel_known_file.close()
                tp_indel_unknown_file.close()
                fp_indel_unknown_file.close()
                
                fp_snp_none_file.close()
                fp_indel_none_file.close()
                
                fn_snp_known_file.close()
                fn_snp_unknown_file.close()
                
                fn_indel_known_file.close()
                fn_indel_unknown_file.close()
                
                fn_snp_known_lowqual_file.close()
                fn_snp_unknown_lowqual_file.close()
                
                fn_indel_known_lowqual_file.close()
                fn_indel_unknown_lowqual_file.close()

                fp_low_qual_file.close()
