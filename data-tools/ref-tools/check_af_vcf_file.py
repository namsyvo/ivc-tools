import sys
import os

vcf_file = open(sys.argv[1])
dir_name, file_name = os.path.split(sys.argv[1])

pos_num = 0
for line in vcf_file:
    if line[0] != '#':
        print
    	print line.split()[7]
        f_gt, s_gt = [], []
        hom_ref, het_01, het_10, hom_var = 0, 0, 0, 0
    	for sample in line.split()[9:]:
            f_gt.append(sample[0])
            s_gt.append(sample[2])
            if sample[0] == '0' and sample[2] == '0':
                hom_ref += 1
            if sample[0] == '0' and sample[2] == '1':
                het_01 += 1
            if sample[0] == '1' and sample[2] == '0':
                het_10 += 1
            if sample[0] == '1' and sample[2] == '1':
                hom_var += 1
        print "#sample, hom_ref, het_01, het_10, hom_var", len(f_gt), hom_ref, het_01, het_10, hom_var

        f_0, f_1, s_0, s_1 = 0, 0, 0, 0
        for f in f_gt:
            if f == '0':
                f_0 += 1
            if f == '1':
                f_1 += 1
        for s in s_gt:
            if s == '0':
                s_0 += 1
            if s == '1':
                s_1 += 1
        print "f_gt_0, f_gt_1, s_gt_0, s_gt_1", f_0, f_1, s_0, s_1

        pos_num += 1
        if pos_num > 5:
            break
