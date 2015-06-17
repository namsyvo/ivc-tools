import sys

t_count, u_count = 0, 0
f = open(sys.argv[1])
for line in f:
    info = line.split("\t")
    if info[1].strip() == "4":
        print info
        u_count += 1
    t_count += 1

print u_count, t_count, float(u_count)/t_count
