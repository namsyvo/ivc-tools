import os
import sys
import json
import re

config_file = open(sys.argv[1])
data = json.load(config_file)
config_file.close()

prog_version = data["ProgVer"]
data_dir = data["DataPath"]["DataDir"]
ref_dir = data["DataPath"]["RefDir"]
genome_fn = data["DataPath"]["GenomeFile"]
result_dir = data["DataPath"]["ResultDir"]
read_fn = data["DataPath"]["ReadPrefixFile"]

confi = float(sys.argv[2])
filter_num = int(sys.argv[3])
result_dn = sys.argv[4]

ref_len = 249250621
ref_para = ['0.70', '0.75', '0.80', '0.85', '0.90', '0.95', '0.96', '0.97', '0.98', '0.99']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
#read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [3, 4, 5, 6, 7, 8, 9, 10]]
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [5]]

ref_path = os.path.join(data_dir, ref_dir)
genome_file = os.path.join(ref_path, genome_fn)

extra_num = 0
intra_num = filter_num - extra_num

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

    trace_info_path = os.path.join(data_dir, result_dir, "alt_mutant_dwgsim", "isc_" + para, result_dn, "snp_num_filter")
    if not os.path.exists(trace_info_path):
        os.makedirs(trace_info_path)
    result_path = os.path.join(data_dir, result_dir, "alt_mutant_dwgsim", "isc_" + para, result_dn)

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                snp, all_snp = {}, {}
                called_snp_file = result_path + "/" + read_fn + "_" + str(rl) + "." + str(err) + "." + str(rn) + ".512.snpcall.32.vcf"
                f = open(called_snp_file)
                for line in f.readlines():
                    value = line.strip().split()
                    if float(value[3]) >= confi:
                        if int(value[0]) - 1 not in snp:
                            snp[int(value[0]) - 1] = []
                        snp[int(value[0]) - 1].append(value)
                    if int(value[0]) - 1 not in all_snp:
                        all_snp[int(value[0]) - 1] = []
                    all_snp[int(value[0]) - 1].append(value)
                print len(snp), len(all_snp)

                sorted_snp = sorted(snp.iteritems(), key=lambda x: float(x[1][0][3]), reverse=True)
                rev_sorted_snp = sorted_snp[::-1]
                file_prefix = trace_info_path + "/" + read_fn + "." + str(rl) + "." + str(err) + "." + str(rn)
                snp_info_file = open(file_prefix + ".snp_info." + str(filter_num) + "." + str(confi) + ".txt", "w")
                snp_info_file.write("rank\tpos\tsnp\tprob\tqual\tnum\taln_dis\tchr_dis\tsnp_indel\ttp_fp\tloc_region\tfilter_region\n")
                rank = 0
                snp_num = 0
                for key, value_arr in rev_sorted_snp:
                    rank += 1
                    for value in value_arr:
                        val_info = "\t".join(value[0:6]) + "\t" + str(abs(int(value[6])))
                        if key in true_snp_comp or key in true_indel_comp:
                            if key in true_snp_comp:
                                snp_num += 1
                                if value[1] == true_snp_comp[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\ttp\tk\textra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\tfp\tk\textra\n")
                            elif key in true_indel_comp:
                                if value[1] == true_indel_comp[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\ttp\tk\textra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\tfp\tk\textra\n")
                        elif key in true_snp_none or key in true_indel_none:
                            if key in true_snp_none:
                                snp_num += 1
                                if value[1] == true_snp_none[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\ttp\tu\textra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\tfp\tu\textra\n")
                            elif key in true_indel_none:
                                if value[1] == true_indel_none[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\ttp\tu\textra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\tfp\tu\textra\n")
                        else:
                            if len(value[1]) == 1:
                                snp_num += 1
                                snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\tfp\to\textra\n")
                            else:
                                snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\tfp\to\textra\n")
                    if snp_num > extra_num:
                        break
                
                intra_snp = rev_sorted_snp[rank:]
                snp_num = 0
                for key, value_arr in intra_snp:
                    rank += 1
                    for value in value_arr:
                        val_info = "\t".join(value[0:6]) + "\t" + str(abs(int(value[6])))
                        if key in true_snp_comp or key in true_indel_comp:
                            if key in true_snp_comp:
                                snp_num += 1
                                if value[1] == true_snp_comp[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\ttp\tk\tintra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\tfp\tk\tintra\n")
                            elif key in true_indel_comp:
                                if value[1] == true_indel_comp[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\ttp\tk\tintra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\tfp\tk\tintra\n")
                        elif key in true_snp_none or key in true_indel_none:
                            if key in true_snp_none:
                                snp_num += 1
                                if value[1] == true_snp_none[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\ttp\tu\tintra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\tfp\tu\tintra\n")
                            elif key in true_indel_none:
                                if value[1] == true_indel_none[key]:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\ttp\tu\tintra\n")
                                else:
                                    snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\tfp\tu\tintra\n")
                        else:
                            if len(value[1]) == 1:
                                snp_num += 1
                                snp_info_file.write(str(rank) + "\t" + val_info + "\tsnp\tfp\to\tintra\n")
                            else:
                                snp_info_file.write(str(rank) + "\t" + val_info + "\tindel\tfp\to\tintra\n")
                    if snp_num > intra_num:
                        break
                
                snp_info_file.close()
