import sys
f=open(sys.argv[1])
snp_count, ins_count, del_count = 0, 0, 0
max_ins_len, max_del_len = 0, 0
ins_len10, ins_len20, ins_len30 = 0, 0, 0
del_len10, del_len20, del_len30 = 0, 0, 0
for line in f:
    if line[0] == '#':
        continue
    tokens = line.strip().split()
    ref, alt = tokens[1], tokens[2]
    if len(ref) == 1 and len(alt) == 1:
        snp_count += 1
    elif len(alt) > 1:
        ins_len = len(alt)
        ins_count += 1
        if max_ins_len < ins_len:
            max_ins_len = ins_len
        if ins_len >= 10:
            ins_len10 += 1
        if ins_len >= 20:
            ins_len20 += 1
        if ins_len >= 30:
            ins_len30 += 1
    elif len(ref) > 1:
        del_len = len(ref)
        del_count += 1
        if max_del_len < del_len:
            max_del_len = del_len
        if del_len >= 10:
            del_len10 += 1
        if del_len >= 20:
            del_len20 += 1
        if del_len >= 30:
            del_len30 += 1
print "snp_count, ins_count, del_count", snp_count, ins_count, del_count
print "ins_len10, ins_len20, ins_len30", ins_len10, ins_len20, ins_len30
print "del_len10, del_len20, del_len30", del_len10, del_len20, del_len30
print "max_ins_len, max_del_len", max_ins_len, max_del_len
