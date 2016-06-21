'''
Exract reads from fastq file with read name obtained from an extracted sam file
'''

import sys

sam_file = open(sys.argv[1])
chr_name = sys.argv[2]

in_fastq_file1 = open("/data/nsvo/test_data/GRCh37_chr1/reads/real_reads/NA12878/ori_reads/ERR194147_sorted_chr1_1.bwa.fastq")
in_fastq_file2 = open("/data/nsvo/test_data/GRCh37_chr1/reads/real_reads/NA12878/ori_reads/ERR194147_sorted_chr1_2.bwa.fastq")

out_fastq_file1 = open(sys.argv[1] + "_chr" + chr_name + "_1.fastq", "w")
out_fastq_file2 = open(sys.argv[1] + "_chr" + chr_name + "_2.fastq", "w")

read_name = []
while True:
    line = sam_file.readline()
    if not line:
        break
    if line[0] != '@':
        read_name.append(line.split("\t")[0])
print read_name

while True:
    line1 = in_fastq_file1.readline()
    line2 = in_fastq_file2.readline()
    if not line1:
        break
    if line1[0] == '@':
        tokens = line1[1:].split("\t")
        rn = tokens[0].split("/")[0]

    if rn in read_name:
        out_fastq_file1.write(line1)
        out_fastq_file1.write(in_fastq_file1.readline())
        out_fastq_file1.write(in_fastq_file1.readline())
        out_fastq_file1.write(in_fastq_file1.readline())

        out_fastq_file2.write(line2)
        out_fastq_file2.write(in_fastq_file2.readline())
        out_fastq_file2.write(in_fastq_file2.readline())
        out_fastq_file2.write(in_fastq_file2.readline())

sam_file.close()
in_fastq_file1.close()
in_fastq_file2.close()
out_fastq_file1.close()
out_fastq_file2.close()
