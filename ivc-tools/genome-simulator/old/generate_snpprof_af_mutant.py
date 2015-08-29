import sys
import os
import random

snp_prof_file = sys.argv[1]
ref_path = sys.argv[2]

snp_prof = {}
for line in open(os.path.join(snp_prof_file)):
    if line[0] != '#':
        info = line.strip().split()
        snp_prof[int(info[1]) - 1] = "\t".join(info[:9])

snp_prof_file = open(os.path.join(ref_path, "snp_prof.vcf"), "w")
with snp_prof_file:
    for key, val in snp_prof.iteritems():
        snp_prof_file.write(val + "\n")


ref_para = ['0.70', '0.75', '0.80', '0.85', '0.90', '0.95', '0.96', '0.97', '0.98', '0.99']

for p in ref_para[0:1]:
    known_variant_file = os.path.join(ref_path, "known_variants.txt")

    snp_prof_comp = {}
    snp_prof_none = {}

    new_var_prob = 1 - float(p)
    for line in open(known_variant_file):
        info = line.strip().split()
        if random.random() < new_var_prob:
            snp_prof_none[int(info[0])] = info[1:]
        else:
            snp_prof_comp[int(info[0])] = info[1:]

    variant_comp_file = open(os.path.join(ref_path, "variant_comp_" + str(p) + ".txt"), "w")
    variant_none_file = open(os.path.join(ref_path, "variant_none_" + str(p) + ".txt"), "w")
    isc_snp_prof_file = open(os.path.join(ref_path, "isc_snp_prof_" + str(p) + ".vcf"), "w")

    for key, val in snp_prof_comp.iteritems():
        variant_comp_file.write(str(key) + "\t")
        for v in val:
            variant_comp_file.write(str(v) + "\t")
        variant_comp_file.write("\n")
        isc_snp_prof_file.write(snp_prof[key] + "\n")

    for key, val in snp_prof_none.iteritems():
        variant_none_file.write(str(key) + "\t")
        for v in val:
            variant_none_file.write(str(v) + "\t")
        variant_none_file.write("\n")

    variant_comp_file.close()
    variant_none_file.close()
    isc_snp_prof_file.close()
