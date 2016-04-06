import sys
f=open(sys.argv[1])
header=f.readline()
data = f.read()
data = data.replace("N", "")
f.close()
f=open(sys.argv[2], "w")
with f:
    f.write(header)
    f.write(data)
f1=open(sys.argv[2])
f2=open(sys.argv[3], "w")
with f2:
    for line in f1:
        if len(line) > 2:
            f2.write(line)
f1.close()
