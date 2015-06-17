import os

prog_path = "/home/nsvo/genome-tools/bwa-0.7.7"
data_path = "/backup/SNP/GRCh37_human"
index_path = data_path + "/bwa-index"

reads_path = "/data/nsvo/Human/SRS006837/reads"
results_path = "/data/nsvo/Human/SRS006837/results"

#reads_path = "/data/nsvo/Human/SRS000204/reads"
#results_path = "/data/nsvo/Human/SRS000204/results"

cmd = "(time " + prog_path + "/bwa bwasw " + index_path + "/GRCh37_chr20.fasta " \
    + reads_path + "/SRR062641_1.fastq " + reads_path + "/SRR062641_2.fastq > " + results_path + "/SRR062641.sam) 2>time-BWA-SRR062641.txt"

#cmd = "time " + prog_path + "/bwa bwasw " + index_path + "/GRCh37_chr1.fasta " \
#    + reads_path + "/SRR352199_1.fastq " + reads_path + "/SRR352199_2.fastq > " + results_path + "/SRR352199-chr1.sam"

print cmd
os.system(cmd)
