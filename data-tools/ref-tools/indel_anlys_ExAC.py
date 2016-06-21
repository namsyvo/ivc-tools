import sys
f=open(sys.argv[1])

snp_count = 0
ins_count = 0
snp_ins_count = 0
snp_del_count = 0
other_del_count = 0
del_count = 0

out1 = open(sys.argv[1] + ".other-del", "w")

for line in f:
    if line[0] == '#':
        continue
    tokens = line.strip().split()
    if tokens[0] != "1":
        continue

    ref, alt = tokens[3], tokens[4]
    allele = alt.split(",")
    other_del = False
    if len(ref) == 1 and len(alt) == 1:
        snp_count += 1
    if len(alt) > 1:
        for a in allele:
            if len(a) > 1:
                ins_count += 1
            else:
                snp_ins_count += 1
    if len(ref) > 1:
        for a in allele:
            if len(a) > 1:
                if ref[0] != a[0] and ref[1:] == a[1:]:
                    snp_del_count += 1
                else:
                    other_del_count += 1
                    other_del = True
            else:
                del_count += 1
    if other_del:
        out1.write(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3] + "\t" + tokens[4] + "\t" + tokens[5] + "\t" + tokens[6] + "\n")
out1.close()

print "snp_count", snp_count
print "ins_count", ins_count
print "snp_ins_count", snp_ins_count
print "snp_del_count", snp_del_count
print "other_del_count", other_del_count
print "del_count", del_count