import sys
import bisect

dbsnp_file = open("/data/nsvo/test_data/GRCh37_chr1/refs/TRIMMED.ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
dbsnp_info = {}
for line in dbsnp_file:
    if line[0] != "#":
        tokens = line.split()
        if tokens[0] == "1":
            dbsnp_info[int(tokens[1]) - 1] = tokens[3:5]
indel = 0
for pos, value in dbsnp_info.iteritems():
    if len(value[0]) > 1 or len(value[1]) > 1:
        indel += 1
        ci = 0
        for p in range(pos-30, pos + 30):
            if p in dbsnp_info and (len(dbsnp_info[p][0]) > 1 or len(dbsnp_info[p][0]) > 1):
                ci += 1
        if ci > 0:
            print ci

print indel
print len(dbsnp_info)

