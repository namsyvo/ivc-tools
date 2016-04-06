import sys

dbsnp_file = open("/backup2/nsvo/variant_calling/Human_data/refs/Whole_Genome/TRIMMED.ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
dbsnp_info, dbsnp_pos = [], []
for line in dbsnp_file:
    if line[0] != "#":
        tokens = line.split()
        if tokens[0] == "1":
            dbsnp_info.append([tokens[1], tokens[3], tokens[4]])
            dbsnp_pos.append(tokens[1])
print len(dbsnp_info)

#var_info_file = open("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NA12878.vcf")
var_info_file = open("/backup2/nsvo/variant_calling/Human_data/refs/NA12878/NISTIntegratedCalls.vcf")
out1 = open(sys.argv[1] + "_1", "w")
out2 = open(sys.argv[1] + "_2", "w")
out3 = open(sys.argv[1] + "_3", "w")
var_info  = []
for line in var_info_file:
    if line[0] != "#":
        tokens = line.split()
        if tokens[0] == "1":
            var_info.append([tokens[1], tokens[3], tokens[4]])
print len(var_info)

c1, c2, c3 = 0, 0, 0
for info in var_info:
    if info[0] in dbsnp_pos:
        if info in dbsnp_info:
            out1.write(str(info) + "\n")
            c1 += 1
        else:
            out2.write(str(info) + "\n")
            c2 += 1
    else:
        out3.write(str(info) + "\n")
        c3 += 1
print c1, c2, c3

out1.close()
out2.close()
out3.close()
