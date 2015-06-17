import sys
fin = open(sys.argv[1])
fref = open(sys.argv[2])
fref.readline()
ref = ''
for line in fref:
	ref += line.strip()
print len(ref)

fout = open(sys.argv[3], 'w')
data = fref.readline
for line in fin:
	if line.strip():
		pos = line.strip().split()[0]
		#if ref[int(pos)] == line.strip().split()[2]:
			#print pos
		if line.strip().split()[1] != line.strip().split()[2]:
			print pos, '\t', line.strip().split()[1], '\t', ref[int(pos)], '\t', '\t'.join(line.strip().split()[2:])
		fout.write(pos + '\t' + ref[int(pos)] + '\n')
fout.close()