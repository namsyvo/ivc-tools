import sys
from random import shuffle

read_file1 = open(sys.argv[1])
read_file2 = open(sys.argv[2])

all_read1, all_read2 = [], []
while True:
    read1, read2 = [], []
    line1 = read_file1.readline()
    line2 = read_file2.readline()
    if not line1 or not line2: break

    read1.append(line1)
    read1.append(read_file1.readline())
    read1.append(read_file1.readline())
    read1.append(read_file1.readline())

    read2.append(line2)
    read2.append(read_file2.readline())
    read2.append(read_file2.readline())
    read2.append(read_file2.readline())

    all_read1.append(read1)
    all_read2.append(read2)

read_file1.close()
read_file2.close()

print len(all_read1), len(all_read2)

id = list(range(len(all_read1)))
shuffle(id)

out_file1 = open(sys.argv[1] + "_shuf", "w")
out_file2 = open(sys.argv[2] + "_shuf", "w")

for i in id:
    out_file1.write(all_read1[i][0])
    out_file1.write(all_read1[i][1])
    out_file1.write(all_read1[i][2])
    out_file1.write(all_read1[i][3])

    out_file2.write(all_read2[i][0])
    out_file2.write(all_read2[i][1])
    out_file2.write(all_read2[i][2])
    out_file2.write(all_read2[i][3])

out_file1.close()
out_file2.close()
