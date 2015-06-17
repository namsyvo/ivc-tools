import sys

snp_index_file = sys.argv[1]

f_in = open(snp_index_file)
f_out = open(snp_index_file + ".add_af", "w")
for line in f_in:
	if line.strip():
		tmp = line.strip().split()
		f_out.write(tmp[0] + "\t")
		l = len(tmp) - 1
		for i in range(1, l + 1):
			f_out.write("%s\t%.5f\t" % (tmp[i], 1.0/l))
		f_out.write("\n")
f_out.close()