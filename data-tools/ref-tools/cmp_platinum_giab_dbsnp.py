import sys
import bisect

bed_file = "/data/nsvo/test_data/GRCh37_chr1/refs/exome_calling_regions.v1.interval_list_chr1.bed"
itv1, itv2 = [], []
with open(bed_file) as f:
    for line in f:
        tokens = line.split("\t")
        itv1.append(int(tokens[1]))
        itv2.append(int(tokens[2]))
print "#intervals", len(itv1)

#var_info_file = open("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NA12878.vcf")
var_info_file = open("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NISTIntegratedCalls.vcf")
var_info  = {}
for line in var_info_file:
    if line[0] != "#":
        tokens = line.split()
        if tokens[0] == "1":
            pos = int(tokens[1]) - 1
            i = bisect.bisect_right(itv1, pos)
            if i: # i.e., itv1[i-1] <= pos
                if pos > itv2[i-1]:
                    continue
            var_info[tokens[1]] = tokens[3:5]
print "#GIAB", len(var_info)

#dbsnp_file = open("/backup2/nsvo/variant_calling/Human_data/refs/Whole_Genome/TRIMMED.ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
dbsnp_file = open("/data/nsvo/test_data/GRCh37_chr1/refs/TRIMMED.ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
#dbsnp_file = open("/data/nsvo/test_data/GRCh37_chr1/refs/ExAC.r0.3.1.sites.vep.vcf_chr1")
dbsnp_info = {}
for line in dbsnp_file:
    if line[0] != "#":
        tokens = line.split()
        if tokens[0] == "1":
            pos = int(tokens[1]) - 1
            i = bisect.bisect_right(itv1, pos)
            if i: # i.e., itv1[i-1] <= pos
                if pos > itv2[i-1]:
                    continue
            dbsnp_info[tokens[1]] = tokens[3:5]
print "#dbsnp", len(dbsnp_info)

out1 = open(sys.argv[1] + "_1", "w")
out2 = open(sys.argv[1] + "_2", "w")
out3 = open(sys.argv[1] + "_3", "w")
c1, c2, c3 = 0, 0, 0
for pos, value in var_info.iteritems():
    if pos in dbsnp_info:
        if value[0] == dbsnp_info[pos][0] and (set(value[1].split(",")).issubset(dbsnp_info[pos][1].split(",")) or set(value[1].split(",")).issuperset(dbsnp_info[pos][1].split(","))):
            out1.write(pos + "\t" + str(value) + "\n")
            c1 += 1
        else:
            out2.write(pos + "\t" + str(value) + "\n")
            c2 += 1
    else:
        out3.write(pos + "\t" + str(value) + "\n")
        c3 += 1
print c1, c2, c3

out1.close()
out2.close()
out3.close()
