import sys
fn = sys.argv[1]
f=open(fn)
len_1 = 0
len_2 = 0
len_3 = 0
len_more = 0
two_snp = 0
three_snp = 0
more_snp = 0
for line in f.readlines():
        tmp = line.strip().split()
        if len(tmp[1:]) == 2 :
                two_snp += 1
        elif len(tmp[1:]) == 3:
                three_snp += 1
        else:
                more_snp += 1
        for snp in tmp[1:]:
                snp_len = len(snp) - 1
       		if snp_len == 1:
                        len_1 += 1
                elif snp_len == 2:
                        len_2 += 1
                elif snp_len == 3:
                        len_3 += 1
                elif snp_len > 3:
                        len_more += 1

print "indel_len = 1: ", len_1
print "indel_len = 2: ", len_2
print "indel_len = 3: ", len_3
print "indel_len > 3: ", len_more
print "total indel: ", len_1 + len_2 + len_3 + len_more

print "two_snp: ", two_snp
print "three_snp: ", three_snp
print "more_snp: ", more_snp
