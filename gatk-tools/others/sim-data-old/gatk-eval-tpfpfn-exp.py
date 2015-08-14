import os
import sys

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
ref_para = ['0.0000', '0.0825', '0.1650', '0.2475', '0.3300']
read_lens = [100]
seq_errs = ['0.00015-0.0015']
read_nums = [cov*ref_len/(2*read_lens[0]) for cov in range(1, 6)]

ref_path = data_path + "/refs"
genome_file = ref_path + "/" + genome_fn

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

    trace_info_path = os.path.join(data_path, "results/sim-reads/", "mutate-" + para + "-dwgsim", "gatk", "trace-info")
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

    result_path = os.path.join(data_path, "results/sim-reads/", "mutate-" + para + "-dwgsim", "gatk")

    for rl in read_lens:
        for err in seq_errs:
            for rn in read_nums[4:5]:                
                snp = {}
                called_snp_file = result_path + "/" + read_fn + "-" + str(rl) + "." + str(err) + "." + str(rn) + ".bwa.vcf"
                f = open(called_snp_file)
                for line in f.readlines():
                    if line[0] != '#':
                        value = line.strip().split()
                        if float(value[5]) >= confi:
                            info = value[7].split(';')
                            for elem in info:
                            	if elem.split('=')[0] == "DP":
                            		dp = elem.split('=')[1]
                            snp[int(value[1]) - 1] = [value[3], value[4], value[5], dp]
                
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

                for key, value in true_snp_comp.iteritems():
                    if key not in snp:
                        fn_snp_comp_file.write(str(key) + "\t" + value + "\n")
                for key, value in true_snp_part.iteritems():
                    if key not in snp:
                        fn_snp_part_file.write(str(key) + "\t" + value + "\n")
                for key, value in true_snp_none.iteritems():
                    if key not in snp:
                        fn_snp_none_file.write(str(key) + "\t" + value + "\n")

                for key, value in true_indel_comp.iteritems():
                    if key not in snp:
                        fn_indel_comp_file.write(str(key) + "\t" + value + "\n")
                for key, value in true_indel_part.iteritems():
                    if key not in snp:
                        fn_indel_part_file.write(str(key) + "\t" + value + "\n")
                for key, value in true_indel_none.iteritems():
                    if key not in snp:
                        fn_indel_none_file.write(str(key) + "\t" + value + "\n")

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
