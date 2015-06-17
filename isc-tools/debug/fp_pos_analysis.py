import sys

f= open(sys.argv[1])
f.readline()
for line in f:
	tmp = line.strip().split()
	pos, chr_diff, header = int(tmp[0]), int(tmp[8]), tmp[15]
	#print pos, chr_diff, header
	tmp = header.split("_")
	pos1, pos2 = int(tmp[2]), int(tmp[3])
	#print pos1, pos2
	diff = 0
	if header[len(header) - 1] == '1':
		diff = abs(pos1 + chr_diff - pos)
	else:
		diff = abs(pos2 + chr_diff - pos)
	print diff