'''
Generate variant profile for IVC as part of set of whole variants
Usage:
    python gen_var_prof_af_sid_mutant.py [var_prof_file_name (e.g., dbSNP)] [ref_path (store mutant variants)]
Output:
    + known_var_(prob_num).txt: set of variants which are known with ISC
    + unknown_var_(prob_num).txt: set of variants which are unknown with ISC
    + isc_var_prof_(prob_num).txt: variant profile which is used as input of ISC (correspoding to known variants)
'''

import sys
import os
import random

if len(sys.argv) == 1:
    var_prof_fn = "/data/nsvo/test-data/GRCh37_chr1/refs/TRIMMED.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.diffcontigname.vcf"
    ref_path = "/data/nsvo/test-data/GRCh37_chr1/refs/af_sid_mutant"
else:
    var_prof_fn = sys.argv[1]
    ref_path = sys.argv[2]

var_prof = {}
for line in open(os.path.join(var_prof_fn)):
    if line[0] != '#':
        info = line.strip().split()
        var_prof[int(info[1]) - 1] = line

known_var_prob_para = ['0.00', '0.50', '0.60', '0.70', '0.80', '0.90', '1.00']

for known_var_prob in known_var_prob_para[6:]:
    known_var_file = open(os.path.join(ref_path, "known_var_" + known_var_prob + ".txt"), "w")
    unknown_var_file = open(os.path.join(ref_path, "unknown_var_" + known_var_prob + ".txt"), "w")
    ivc_var_prof_file = open(os.path.join(ref_path, "ivc_var_prof_" + known_var_prob + ".vcf"), "w")

    mut_var_file = os.path.join(ref_path, "mut_var.txt")
    for line in open(mut_var_file):
        if line[0] == '#':
            continue
        info = line.strip().split()
        if random.random() < float(known_var_prob):
            known_var_file.write(line)
            ivc_var_prof_file.write(var_prof[int(info[0])])
        else:
            unknown_var_file.write(line)
    known_var_file.close()
    unknown_var_file.close()
    ivc_var_prof_file.close()
