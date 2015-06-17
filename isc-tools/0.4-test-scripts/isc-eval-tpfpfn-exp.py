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

ref_len = 249250621
ref_para = ['0.0000', '0.0825', '0.1650', '0.2475', '0.3300']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in [25]]

ref_path = os.path.join(data_path, ref_dir)
genome_file = os.path.join(ref_path, genome_fn)

result_dn = sys.argv[3]

for para in ref_para[4:5]:

    true_snp_comp, true_indel_comp, true_snp_part, true_indel_part, true_snp_none, true_indel_none = {}, {}, {}, {}, {}, {}

    variant_comp_file = os.path.join(ref_path, "mutate-" + para, "variant_comp.txt")
    variant_part_file = os.path.join(ref_path, "mutate-" + para, "variant_part.txt")
    variant_none_file = os.path.join(ref_path, "mutate-" + para, "variant_none.txt")

    with open(variant_comp_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                if len(value[1]) == 1 and value[1] != ".":
                    true_snp_comp[int(value[0])] = value[1]
                else:
                    true_indel_comp[int(value[0])] = value[1]

    with open(variant_part_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                if len(value[1]) == 1 and value[1] != ".":
                    true_snp_part[int(value[0])] = value[1]
                else:
                    true_indel_part[int(value[0])] = value[1]

    with open(variant_none_file) as f:
        for line in f.readlines():
            if line.strip():
                value = line.strip().split()
                if len(value[1]) == 1 and value[1] != ".":
                    true_snp_none[int(value[0])] = value[1]
                else:
                    true_indel_none[int(value[0])] = value[1]

    trace_info_path = os.path.join(data_path, result_dir, "mutate-" + para + "-dwgsim", "isc", "debug", "trace-info", result_dn)
    if not os.path.exists(trace_info_path):
        os.makedirs(trace_info_path)
    file_prefix = trace_info_path + "/" + read_fn + "-" + str(read_lens[0]) + "." + str(seq_errs[0])

    tp_snp_none_file = open(file_prefix + ".tp_snp_none." + str(confi) + ".txt", "w")
    fp_snp_none_file = open(file_prefix + ".fp_snp_none." + str(confi) + ".txt", "w")
    tp_indel_none_file = open(file_prefix + ".tp_indel_none." + str(confi) + ".txt", "w")
    fp_indel_none_file = open(file_prefix + ".fp_indel_none." + str(confi) + ".txt", "w")

    tp_snp_part_file = open(file_prefix + ".tp_snp_part." + str(confi) + ".txt", "w")
    fp_snp_part_file = open(file_prefix + ".fp_snp_part." + str(confi) + ".txt", "w")
    tp_indel_part_file = open(file_prefix + ".tp_indel_part." + str(confi) + ".txt", "w")
    fp_indel_part_file = open(file_prefix + ".fp_indel_part." + str(confi) + ".txt", "w")

    tp_snp_comp_file = open(file_prefix + ".tp_snp_comp." + str(confi) + ".txt", "w")
    fp_snp_comp_file = open(file_prefix + ".fp_snp_comp." + str(confi) + ".txt", "w")
    tp_indel_comp_file = open(file_prefix + ".tp_indel_comp." + str(confi) + ".txt", "w")
    fp_indel_comp_file = open(file_prefix + ".fp_indel_comp." + str(confi) + ".txt", "w")

    fp_snp_other_file = open(file_prefix + ".fp_snp_other." + str(confi) + ".txt", "w")
    fp_indel_other_file = open(file_prefix + ".fp_indel_other." + str(confi) + ".txt", "w")

    fn_snp_comp_file = open(file_prefix + ".fn_snp_comp." + str(confi) + ".txt", "w")
    fn_snp_part_file = open(file_prefix + ".fn_snp_part." + str(confi) + ".txt", "w")
    fn_snp_none_file = open(file_prefix + ".fn_snp_none." + str(confi) + ".txt", "w")

    fn_indel_comp_file = open(file_prefix + ".fn_indel_comp." + str(confi) + ".txt", "w")
    fn_indel_part_file = open(file_prefix + ".fn_indel_part." + str(confi) + ".txt", "w")
    fn_indel_none_file = open(file_prefix + ".fn_indel_none." + str(confi) + ".txt", "w")

    fn_snp_comp_callgatk_file = open(file_prefix + ".fn_snp_comp_callgatk." + str(confi) + ".txt", "w")
    fn_snp_part_callgatk_file = open(file_prefix + ".fn_snp_part_callgatk." + str(confi) + ".txt", "w")
    fn_snp_none_callgatk_file = open(file_prefix + ".fn_snp_none_callgatk." + str(confi) + ".txt", "w")

    fn_indel_comp_callgatk_file = open(file_prefix + ".fn_indel_comp_callgatk." + str(confi) + ".txt", "w")
    fn_indel_part_callgatk_file = open(file_prefix + ".fn_indel_part_callgatk." + str(confi) + ".txt", "w")
    fn_indel_none_callgatk_file = open(file_prefix + ".fn_indel_none_callgatk." + str(confi) + ".txt", "w")

    fp_all_file = open(file_prefix + ".fp_all." + str(confi) + ".txt", "w")

    result_path = os.path.join(data_path, result_dir, "mutate-" + para + "-dwgsim", "isc", "debug", result_dn)
    gatk_result_path = os.path.join(data_path, result_dir, "mutate-" + para + "-dwgsim", "gatk")

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums:
                snp = {}
                all_snp = {}
                called_snp_file = result_path + "/" + read_fn + "-" + str(rl) + "." + str(err) + "." + str(rn) + ".4096.snpcall.32.vcf"
                f = open(called_snp_file)
                for line in f.readlines():
                    value = line.strip().split()
                    if float(value[2]) >= confi:
                        snp[int(value[0]) - 1] = value
                    all_snp[int(value[0]) - 1] = value
                
                gatk_snp = {}
                gatk_called_snp_file = gatk_result_path + "/" + read_fn + "-" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.vcf"
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

                for key, value in snp.iteritems():
                    val_info = str(key) + "\t" + "\t".join(value[1:]) + "\n"
                    if key in true_snp_comp or key in true_indel_comp:
                        if key in true_snp_comp:
                            if value[1] == true_snp_comp[key]:
                                tp_snp_comp_file.write(val_info)
                            else:
                                fp_snp_comp_file.write(val_info)
                        elif key in true_indel_comp:
                            if value[1] == true_indel_comp[key]:
                                tp_indel_comp_file.write(val_info)
                            else:
                                fp_indel_comp_file.write(val_info)
                    elif key in true_snp_part or key in true_indel_part:
                        if key in true_snp_part:
                            if value[1] == true_snp_part[key]:
                                tp_snp_part_file.write(val_info)
                            else:
                                fp_snp_part_file.write(val_info)
                        elif key in true_indel_part:
                            if value[1] == true_indel_part[key]:
                                tp_indel_part_file.write(val_info)
                            else:
                                fp_indel_part_file.write(val_info)
                    elif key in true_snp_none or key in true_indel_none:
                        if key in true_snp_none:
                            if value[1] == true_snp_none[key]:
                                tp_snp_none_file.write(val_info)
                            else:
                                fp_snp_none_file.write(val_info)
                        elif key in true_indel_none:
                            if value[1] == true_indel_none[key]:
                                tp_indel_none_file.write(val_info)
                            else:
                                fp_indel_none_file.write(val_info)
                    else:
                            if len(value[1]) == 1:
                                fp_snp_other_file.write(val_info)
                            else:
                                fp_indel_other_file.write(val_info)

                #Print FN info
                for key, value in true_snp_comp.iteritems():
                    if key not in snp:
                        fn_snp_comp_file.write(str(key) + "\t" + value + "\t")
                        if key in all_snp:
                            fn_snp_comp_file.write("\t".join(all_snp[key]))
                            if key in gatk_snp:
                                fn_snp_comp_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")
                        fn_snp_comp_file.write("\n")

                for key, value in true_snp_part.iteritems():
                    if key not in snp:
                        fn_snp_part_file.write(str(key) + "\t" + value + "\t")
                        if key in all_snp:
                            fn_snp_part_file.write("\t".join(all_snp[key]))
                            if key in gatk_snp:
                                fn_snp_part_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")
                        fn_snp_part_file.write("\n")

                for key, value in true_snp_none.iteritems():
                    if key not in snp:
                        fn_snp_none_file.write(str(key) + "\t" + value + "\t")
                        if key in all_snp:
                            fn_snp_none_file.write("\t".join(all_snp[key]))
                            if key in gatk_snp:
                                fn_snp_none_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")
                        fn_snp_none_file.write("\n")

                for key, value in true_indel_comp.iteritems():
                    if key not in snp:
                        fn_indel_comp_file.write(str(key) + "\t" + value + "\t")
                        if key in all_snp:
                            fn_indel_comp_file.write("\t".join(all_snp[key]))
                            if key in gatk_snp:
                                fn_indel_comp_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")
                        fn_indel_comp_file.write("\n")

                for key, value in true_indel_part.iteritems():
                    if key not in snp:
                        fn_indel_part_file.write(str(key) + "\t" + value + "\t")
                        if key in all_snp:
                            fn_indel_part_file.write("\t".join(all_snp[key]))
                            if key in gatk_snp:
                                fn_indel_part_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")
                        fn_indel_part_file.write("\n")

                for key, value in true_indel_none.iteritems():
                    if key not in snp:
                        fn_indel_none_file.write(str(key) + "\t" + value + "\t")
                        if key in all_snp:
                            fn_indel_none_file.write("\t".join(all_snp[key]))
                            if key in gatk_snp:
                                fn_indel_none_callgatk_file.write(str(key) + "\t" + "\t".join(gatk_snp[key]) + "\n")
                        fn_indel_none_file.write("\n")

                for key, value in all_snp.iteritems():
                    if key not in true_snp_comp and key not in true_snp_part and key not in true_snp_none \
                        and key not in true_indel_comp and key not in true_indel_part and key not in true_indel_none:
                        fp_all_file.write(str(key) + "\t" + "\t".join(value) + "\n")
                        if key in gatk_snp:
                            fp_all_file.write(str(key) + "\t".join(gatk_snp[key]) + "\n")

tp_snp_comp_file.close()
fp_snp_comp_file.close()
tp_snp_part_file.close()
fp_snp_part_file.close()
tp_snp_none_file.close()
fp_snp_none_file.close()

tp_indel_comp_file.close()
fp_indel_comp_file.close()
tp_indel_part_file.close()
fp_indel_part_file.close()
tp_indel_none_file.close()
fp_indel_none_file.close()

fp_snp_other_file.close()
fp_indel_other_file.close()

fn_snp_comp_file.close()
fn_snp_part_file.close()
fn_snp_none_file.close()

fn_indel_comp_file.close()
fn_indel_part_file.close()
fn_indel_none_file.close()

fn_snp_comp_callgatk_file.close()
fn_snp_part_callgatk_file.close()
fn_snp_none_callgatk_file.close()

fn_indel_comp_callgatk_file.close()
fn_indel_part_callgatk_file.close()
fn_indel_none_callgatk_file.close()

fp_all_file.close()
