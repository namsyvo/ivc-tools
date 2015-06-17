import sys
import numpy
ins_size_arr = []
f = open(sys.argv[1])
f2 = open(sys.argv[1] + ".isize.txt", "w")
rev_count = 0
for line in f:
    if line[0:3] == "@gi":
        info = line.split('_')
        ins_size = int(info[2]) - int(info[3])
        ins_size_arr.append(abs(ins_size))
        f2.write(str(abs(ins_size)) + "\n")
        if ins_size < 0:
            rev_count += 1
f.close()
f2.close()
print "Read num", len(ins_size_arr)
print "Rev count", rev_count
a = numpy.array(ins_size_arr)
print "Mean of Ins size:", numpy.mean(a)
print "Std Dev of Ins size", numpy.std(a)
