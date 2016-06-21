import sys
fin = open(sys.argv[1])
#fout = open(sys.argv[1] + "_chr1.bed", "w")
fout = open(sys.argv[1] + "_chr1", "w")
chr = sys.argv[2]
for line in fin:
    if line.split("\t")[0] == chr:
        fout.write(line)
fin.close()
fout.close()

