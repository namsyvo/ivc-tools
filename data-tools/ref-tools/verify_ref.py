'''
ref_path = "/home/SNP/NC_016088.1_soybean/index"
ref = []
with open(ref_path + "/genomestar.txt") as f:
	line = f.readline()
	ref.extend(line.strip())
print len(ref)
print ref[0], ref[len(ref)-1]

N_count = 0
for c in ref:
	if c == 'N':
		N_count += 1
print N_count
'''
ref_path = "/home/SNP/GRCh37_human/index/"
ref = []
with open(ref_path + "/genomestar.txt") as f:
	line = f.readline()
	ref.extend(line.strip())
print len(ref)
print ref[0], ref[len(ref)-1]

N_count = 0
for c in ref:
	if c == 'N':
		N_count += 1
print N_count
