import sys
import random

snp_file = sys.argv[1]
contigname = sys.argv[2]
strand = sys.argv[3]

inf = open(snp_file)
outf = open(snp_file + ".dwgsim", "w")
d1, d2 = 0, 0
for line in inf:
	if line[0] != '#':
		new_line = []
		info = line.split()
		new_line.append(contigname)
		new_line.extend(info[1:4])
                if info[4] == ".":
                        d += 1
                if info[4] == "<DEL>":
                        d2 += 1
                        continue
		#if random.random() < 1:
		new_line.append(info[4])
		#else:
		#	new_line.append(info[3])
		new_line.extend(info[5:7])
		#for future: AF in info[7] for new_line need to be modified based on the selected mutation
		new_line.append("pl=" + strand + ";" + info[7])
		new_line.append(info[8])
		outf.write("\t".join(new_line) + "\n")
	elif line[0:2] != "##":
		info = line.split()
		outf.write("\t".join(info[:9]) + "\n")
	else:
		outf.write(line)
print d1, d2
outf.close()
inf.close()
