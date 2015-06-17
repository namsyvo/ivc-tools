import os
import sys
import json

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

prog_version = data["ProgVer"]
data_path = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
genome_fn = data["DataPath"]["GenomeFile"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]

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

ref_path = os.path.join(data_path, ref_dir)
genome_file = os.path.join(ref_path, genome_fn)

result_dn = sys.argv[4]

for para in ref_para:

    true_snp_comp, true_indel_comp, true_snp_none, true_indel_none, var_prof = {}, {}, {}, {}, {}

    variant_comp_file = os.path.join(ref_path, "variant_comp_" + para + ".txt")
    variant_none_file = os.path.join(ref_path, "variant_none_" + para + ".txt")
    variant_prof_file = os.path.join(ref_path, "snp_prof.vcf")

    with open(variant_prof_file) as f:
        for line in f.readlines():
            if line.strip() and line[0] != "#":
                value = line.strip().split()
                var_prof[int(value[1]) - 1] = value[3:5]

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

    fpfntp_info_path = os.path.join(data_path, result_dir, "isc_" + para, result_dn, "fpfntp_info")
    if not os.path.exists(fpfntp_info_path):
        os.makedirs(fpfntp_info_path)

    result_path = os.path.join(data_path, result_dir, "isc_" + para, result_dn)
    gatk_result_path = os.path.join(data_path, result_dir, "gatk")

    header = "pos\tcall_snp\tqual\tprob\tbase_num\tbase_qual\tchr_dis\tchr_diff\t-log(aln_prob)\t-log(pair_prob)\ts_pos1\tbranch1\ts_pos2\tbranch2\tread_header\taln_base\taln_base_num\n"
    fp_header = "pos\ttrue_snp\tcall_snp\tqual\tprob\tbase_num\tbase_qual\tchr_dis\tchr_diff\t-log(aln_prob)\t-log(pair_prob)\ts_pos1\tbranch1\ts_pos2\tbranch2\tread_header\taln_base\taln_base_num\n"

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                fn_part = read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn)
                gatk_snp = {}
                gatk_called_snp_file = gatk_result_path + "/" + fn_part + ".bwa.vcf"
                f = open(gatk_called_snp_file)
                for line in f.readlines():
                    if line[0] != '#':
                        value = line.strip().split()
                        if float(value[5]) >= 20.0:
                            info = value[7].split(';')
                            for elem in info:
                                if elem.split('=')[0] == "DP":
                                    dp = elem.split('=')[1]
                            gatk_snp[int(value[1]) - 1] = [value[3], value[4], value[5], dp]
                f.close()
                print "GATK", len(gatk_snp)

                file_prefix = fpfntp_info_path + "/" + fn_part

                tp_snp_none_file = open(file_prefix + ".tp_snp_none." + str(confi) + ".txt", "w")
                tp_snp_none_file.write(header)
                fp_snp_none_file = open(file_prefix + ".fp_snp_none." + str(confi) + ".txt", "w")
                fp_snp_none_file.write(fp_header)

                tp_indel_none_file = open(file_prefix + ".tp_indel_none." + str(confi) + ".txt", "w")
                tp_indel_none_file.write(header)
                fp_indel_none_file = open(file_prefix + ".fp_indel_none." + str(confi) + ".txt", "w")
                fp_indel_none_file.write(fp_header)

                tp_snp_comp_file = open(file_prefix + ".tp_snp_comp." + str(confi) + ".txt", "w")
                tp_snp_comp_file.write(header)
                fp_snp_comp_file = open(file_prefix + ".fp_snp_comp." + str(confi) + ".txt", "w")
                fp_snp_comp_file.write(fp_header)
                tp_indel_comp_file = open(file_prefix + ".tp_indel_comp." + str(confi) + ".txt", "w")
                tp_indel_comp_file.write(header)
                fp_indel_comp_file = open(file_prefix + ".fp_indel_comp." + str(confi) + ".txt", "w")
                fp_indel_comp_file.write(fp_header)

                fp_snp_other_file = open(file_prefix + ".fp_snp_other." + str(confi) + ".txt", "w")
                fp_snp_other_file.write(header)
                fp_indel_other_file = open(file_prefix + ".fp_indel_other." + str(confi) + ".txt", "w")
                fp_indel_other_file.write(header)

                fn_snp_comp_file = open(file_prefix + ".fn_snp_comp." + str(confi) + ".txt", "w")
                fn_snp_comp_file.write("pos\tsnp\n")
                fn_snp_none_file = open(file_prefix + ".fn_snp_none." + str(confi) + ".txt", "w")
                fn_snp_none_file.write("pos\tsnp\n")

                fn_indel_comp_file = open(file_prefix + ".fn_indel_comp." + str(confi) + ".txt", "w")
                fn_indel_comp_file.write("pos\tsnp\n")
                fn_indel_none_file = open(file_prefix + ".fn_indel_none." + str(confi) + ".txt", "w")
                fn_indel_none_file.write("pos\tsnp\n")

                fn_snp_comp_lowqual_file = open(file_prefix + ".fn_snp_comp_lowqual." + str(confi) + ".txt", "w")
                fn_snp_none_lowqual_file = open(file_prefix + ".fn_snp_none_lowqual." + str(confi) + ".txt", "w")

                fn_indel_comp_lowqual_file = open(file_prefix + ".fn_indel_comp_lowqual." + str(confi) + ".txt", "w")
                fn_indel_none_lowqual_file = open(file_prefix + ".fn_indel_none_lowqual." + str(confi) + ".txt", "w")

                fn_snp_comp_callgatk_file = open(file_prefix + ".fn_snp_comp_callgatk." + str(confi) + ".txt", "w")
                fn_snp_none_callgatk_file = open(file_prefix + ".fn_snp_none_callgatk." + str(confi) + ".txt", "w")

                fn_indel_comp_callgatk_file = open(file_prefix + ".fn_indel_comp_callgatk." + str(confi) + ".txt", "w")
                fn_indel_none_callgatk_file = open(file_prefix + ".fn_indel_none_callgatk." + str(confi) + ".txt", "w")

                fp_low_qual_file = open(file_prefix + ".fp_low_qual." + str(confi) + ".txt", "w")
                fp_low_qual_gatk_file = open(file_prefix + ".fp_low_qual_gatk." + str(confi) + ".txt", "w")

                called_snp, low_qual_snp = {}, {}
                called_snp_file = result_path + "/" + fn_part + ".512.snpcall.32.vcf"
                f = open(called_snp_file)
                for line in f.readlines():
                    value = line.strip().split()
                    key = int(value[0]) - 1
                    if float(value[2]) >= confi:
                        called_snp[key] = value
                        val_info = str(key) + "\t" + "\t".join(value[1:]) + "\n"
                        if key in true_snp_comp or key in true_indel_comp:
                            if value[1] != var_prof[key][0]:
                                if key in true_snp_comp:
                                    if value[1] == true_snp_comp[key]:
                                        tp_snp_comp_file.write(val_info)
                                    else:
                                        fp_snp_comp_file.write(str(key) + "\t" + true_snp_comp[key] + "\t" + "\t".join(value[1:]) + "\n")
                                elif key in true_indel_comp:
                                    if value[1] == true_indel_comp[key]:
                                        tp_indel_comp_file.write(val_info)
                                    else:
                                        fp_indel_comp_file.write(str(key) + "\t" + true_indel_comp[key] + "\t" + "\r".join(value[1:]) + "\n")
                        elif key in true_snp_none or key in true_indel_none:
                            if key in true_snp_none:
                                if value[1] == true_snp_none[key]:
                                    tp_snp_none_file.write(val_info)
                                else:
                                    fp_snp_none_file.write(str(key) + "\t" + true_snp_none[key] + "\t" + "\t".join(value[1:]) + "\n")
                            elif key in true_indel_none:
                                if value[1] == true_indel_none[key]:
                                    tp_indel_none_file.write(val_info)
                                else:
                                    fp_indel_none_file.write(str(key) + "\t" + true_indel_none[key] + "\t" + "\t".join(value[1:]) + "\n")
                        else:
                            if len(value[1]) == 1:
                                fp_snp_other_file.write(val_info)
                            else:
                                fp_indel_other_file.write(val_info)
                    else:
                        low_qual_snp[key] = value

                f.close()
                print "ISC", len(called_snp)

                #Print FN info
                for key, value in true_snp_comp.iteritems():
                    if key not in called_snp and value != var_prof[key][0]:
                        fn_snp_comp_file.write(str(key) + "\t" + value + "\n")
                        if key in low_qual_snp:
                            fn_snp_comp_lowqual_file.write(str(key) + "\t" + "\t".join(low_qual_snp[key]) + "\n")
                        if key in gatk_snp:
                            fn_snp_comp_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")

                for key, value in true_snp_none.iteritems():
                    if key not in called_snp:
                        fn_snp_none_file.write(str(key) + "\t" + value + "\n")
                        if key in low_qual_snp:
                            fn_snp_none_lowqual_file.write(str(key) + "\t" + "\t".join(low_qual_snp[key]) + "\n")
                        if key in gatk_snp:
                            fn_snp_none_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")

                for key, value in true_indel_comp.iteritems():
                    if key not in called_snp and value != var_prof[key][0]:
                        fn_indel_comp_file.write(str(key) + "\t" + value + "\n")
                        if key in low_qual_snp:
                            fn_indel_comp_lowqual_file.write(str(key) + "\t" + "\t".join(low_qual_snp[key]) + "\n")
                        if key in gatk_snp:
                            fn_indel_comp_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")

                for key, value in true_indel_none.iteritems():
                    if key not in called_snp:
                        fn_indel_none_file.write(str(key) + "\t" + value + "\n")
                        if key in low_qual_snp:
                            fn_indel_none_lowqual_file.write(str(key) + "\t" + "\t".join(low_qual_snp[key]) + "\n")
                        if key in gatk_snp:
                            fn_indel_none_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")

                for key, value in low_qual_snp.iteritems():
                    if key not in true_snp_comp and key not in true_snp_none \
                        and key not in true_indel_comp and key not in true_indel_none and value != var_prof[key][0]:
                        fp_low_qual_file.write(str(key) + "\t" + "\t".join(value) + "\n")
                        if key in gatk_snp:
                            fp_low_qual_gatk_file.write(str(key) + "\t".join(gatk_snp[key]) + "\n")

                tp_snp_comp_file.close()
                fp_snp_comp_file.close()
                tp_snp_none_file.close()
                fp_snp_none_file.close()
                
                tp_indel_comp_file.close()
                fp_indel_comp_file.close()
                tp_indel_none_file.close()
                fp_indel_none_file.close()
                
                fp_snp_other_file.close()
                fp_indel_other_file.close()
                
                fn_snp_comp_file.close()
                fn_snp_none_file.close()
                
                fn_indel_comp_file.close()
                fn_indel_none_file.close()
                
                fn_snp_comp_lowqual_file.close()
                fn_snp_none_lowqual_file.close()
                
                fn_indel_comp_lowqual_file.close()
                fn_indel_none_lowqual_file.close()

                fn_snp_comp_callgatk_file.close()
                fn_snp_none_callgatk_file.close()
                
                fn_indel_comp_callgatk_file.close()
                fn_indel_none_callgatk_file.close()

                fp_low_qual_file.close()
                fp_low_qual_gatk_file.close()
