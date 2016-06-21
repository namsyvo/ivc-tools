import sys

sam_file = open(sys.argv[1])
chr_name = sys.argv[2]

fastq_file1 = open(sys.argv[1] + "_chr" + chr_name + "_1.fastq", "w")
fastq_file2 = open(sys.argv[1] + "_chr" + chr_name + "_2.fastq", "w")

'''
fastq_file1 = open("/data/nsvo/test_data/GRCh37_chr22/reads/real_reads/NA12878/ERR194147_sorted_chr22_1.bwa.fastq", "w")
fastq_file2 = open("/data/nsvo/test_data/GRCh37_chr22/reads/real_reads/NA12878/ERR194147_sorted_chr22_2.bwa.fastq", "w")
'''

prev_pos = sam_file.tell()
while True:
    line = sam_file.readline()
    if not line:
        break
    tokens1 = line.split()
    prev_pos = sam_file.tell()
    line = sam_file.readline()
    if not line:
        break
    tokens2 = line.split()

    if tokens1[0] == tokens2[0]:
        if tokens1[2] == tokens2[2]:
            if tokens1[2] == chr_name:
                fastq_file1.write("@" + tokens1[0] + "\n")
                fastq_file1.write(tokens1[9] + "\n")
                fastq_file1.write("+\n")
                fastq_file1.write(tokens1[10] + "\n")
                fastq_file2.write("@" + tokens2[0] + "\n")
                fastq_file2.write(tokens2[9] + "\n")
                fastq_file2.write("+\n")
                fastq_file2.write(tokens2[10] + "\n")
        else:
            print "Not same chr", tokens1[2], tokens2[2]
    else:
        print "Not paired end", tokens1[0], tokens2[0]
        sam_file.seek(prev_pos)

sam_file.close()
fastq_file1.close()
fastq_file2.close()
