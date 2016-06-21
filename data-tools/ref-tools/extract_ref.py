import sys
d = []
f = open(sys.argv[1])
data = f.readlines()
d = ""
for line in data:
    if line[0:2] == '>2':
        break
    else:
        if line[0] == '>':
            continue
        d += line.strip()
print len(d)
pos = int(sys.argv[2])
print str(d[pos:pos+100])
