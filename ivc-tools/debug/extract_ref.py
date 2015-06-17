import sys
fref = open(sys.argv[1])
ref = ''
for line in fref:
	if line[0] != '>':
		ref += line.strip()
s_pos = int(sys.argv[2])
dis = int(sys.argv[3])
read = ref[s_pos : s_pos + dis]
print read
