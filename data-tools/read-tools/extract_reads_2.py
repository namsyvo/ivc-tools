import sys
import random

for cov in [1, 2, 5, 10, 20]:

    in_file1 = open("/data/nsvo/test_data/GRCh37_chr1/reads/real_reads/NA12878/ERR194147_sorted_chr1_1.bwa.fastq")
    in_file2 = open("/data/nsvo/test_data/GRCh37_chr1/reads/real_reads/NA12878/ERR194147_sorted_chr1_2.bwa.fastq")

    out_file1 = open("/data/nsvo/test_data/GRCh37_chr1/reads/real_reads/NA12878/ERR194147_sorted_" + str(cov) + "x_chr1_1.bwa.fastq", "w")
    out_file2 = open("/data/nsvo/test_data/GRCh37_chr1/reads/real_reads/NA12878/ERR194147_sorted_" + str(cov) + "x_chr1_2.bwa.fastq", "w")

    while True:
        line1 = in_file1.readline()
        line2 = in_file2.readline()
        if not line1 or not line2:
            break
        if random.random() < 0.02 * cov:
            out_file1.write(line1)
            out_file2.write(line2)
            out_file1.write(in_file1.readline())
            out_file2.write(in_file2.readline())
            out_file1.write(in_file1.readline())
            out_file2.write(in_file2.readline())
            out_file1.write(in_file1.readline())
            out_file2.write(in_file2.readline())
        else:
            in_file1.readline()
            in_file2.readline()
            in_file1.readline()
            in_file2.readline()
            in_file1.readline()
            in_file2.readline()

    out_file1.close()
    out_file2.close()

    in_file1.close()
    in_file2.close()
