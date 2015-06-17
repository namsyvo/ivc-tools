'''
snp_prof_file = open("/data/nsvo/test-data/GRCh37_chr1/refs/af_mutant/isc_snp_prof_0.70.vcf")
known_true_snp_file = open("/data/nsvo/test-data/GRCh37_chr1/refs/af_mutant/variant_comp_0.70.txt")
isc_snp_prof_file = open("/data/nsvo/test-data/GRCh37_chr1/indexes/af_mutant_index/index_0.70/isc_snp_prof_0.70.vcf.idx")

snp_prof, known_true_snp, isc_snp_prof = {}, {}, {}
for line in snp_prof_file:
	info = line.split()
	snp_prof[int(info[1]) - 1] = info[3:5]

for line in isc_snp_prof_file:
	info = line.split()
	isc_snp_prof[int(info[0])] = info[1:]

for line in known_true_snp_file:
	info = line.split()
	known_true_snp[int(info[0])] = info[1]

a, b = 0, 0
for key, value in isc_snp_prof.iteritems():
	a += 1
	if key not in known_true_snp:
		b += 1
		#print key, value

print a, b
'''

true_snp_file = open("/data/nsvo/test-data/GRCh37_chr1/refs/af_mutant/known_variants.txt")

true_snp = {}
for line in true_snp_file:
	info = line.split()
	true_snp[int(info[0])] = info[1]

sorted_snp = sorted(true_snp.iteritems(), key=lambda x: x[0], reverse=False)

for key, value in sorted_snp:
	if len(value) > 1:
		print key, value
	if key > 1070000:
		break