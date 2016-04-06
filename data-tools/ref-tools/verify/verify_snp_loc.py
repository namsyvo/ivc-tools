import sys
fn, l = sys.argv[1], int(sys.argv[2])
f=open(fn)
n_indel = 0
del_num = 0
ins_num = 0
long_indel_num = 0
two_alle_snp = 0
total_snp = 0
for line in f.readlines():
	tmp = line.strip().split()
	flag_indel = False
	max_len = len(tmp[1])
	if len(tmp[1]) > 1:
		del_num += 1
	if len(tmp[2]) > 1:
		ins_num += 1
	if len(tmp[1]) > 1 and len(tmp[2]) > 1 :
		long_indel_num += 1

	if len(tmp[1:]) > 2:
		print tmp[1:]
	for snp in tmp[1:]:
		if len(snp) > 1:
			flag_indel = True
		if max_len < len(snp):
			max_len = len(snp)
	if flag_indel:
		n_indel += 1
		#print max_len
	else:
		total_snp += 1
		if len(tmp[1:]) == 2:
			two_alle_snp += 1		

print "indel num:", n_indel
print "del num:", del_num
print "ins num:", ins_num
print "long indel num:", long_indel_num
print "snp with 2 alleles:", two_alle_snp
print "total snp:", total_snp