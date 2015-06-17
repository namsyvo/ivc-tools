import os

data_path = "/backup/SNP/GRCh37_human"
index_path = data_path + "/bwa-index"

reads_path = "/data/nsvo/Human/SRS000204/reads"
results_path = "/data/nsvo/Human/SRS000204/results"

prog_path = "/home/nsvo/genome-tools/bwa-0.7.7"

cmd = "time " + prog_path + "/bwa bwasw " + index_path + "/GRCh37_chr20.fasta " \
    + reads_path + "/SRR352199_1.fastq " + reads_path + "/SRR352199_2.fastq > " + results_path + "/SRR352199-chr20.sam"
print cmd
os.system(cmd)

ref_path = "/backup/SNP/GRCh37_human/genomes"
sam_path = "/data/nsvo/Human/SRS000204/results"
tool_path = "/home/nsvo/genome-tools"
script_path = "/home/nsvo/workspace/isc-scripts/gatktools/real-data"

cmd = "time " + script_path + "/gatk-preprocess-real-data.sh " + ref_path + "/GRCh37_chr20 " \
	+ sam_path + "/SRR352199-chr20 " + tool_path
os.system(cmd)

ref_path = "/backup/SNP/GRCh37_human/genomes"
sam_path = "/data/nsvo/Human/SRS000204/results"
results_path = "/data/nsvo/Human/SRS000204/results"
script_path = "/home/nsvo/workspace/isc-scripts/gatktools/real-data"
tool_path = "/home/nsvo/genome-tools"

cmd = "time " + script_path + "/gatk-callsnp.sh " + ref_path + "/GRCh37_chr20.fasta " \
	+ sam_path + "/SRR352199-chr20_sorted_RG.bam " \
	+ ref_path + "/ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes-diffcontigname.vcf " \
	+ results_path + "/SRR352199-chr20.vcf " + tool_path
os.system(cmd)

result_path = "/data/nsvo/Human/SRS000204/results"
script_path = "/home/nsvo/workspace/isc-scripts/gatktools"

cmd = "go run " + script_path + "/extract-alt-vcf.go " + result_path + "/SRR352199-chr20.vcf > " \
	+ result_path + "/SRR352199-chr20.vcf.txt"
os.system(cmd)

import sys

db_snp = {}
db_indel = {}
snp_file = "/backup/SNP/GRCh37_human/index/SNPLocation.txt"
with open(snp_file) as f:
        for line in f.readlines():
                if line.strip():
                        value=line.strip().split()
                        db_snp[int(value[0])] = value[1:]
                        if len(value[1]) > 1:
                                db_indel[int(value[0])] = value[1:]

result_path = "/data/nsvo/Human/SRS000204/results"

snp = {}
called_snp_file = result_path + "/SRR352199-chr20.vcf.txt"
with open(called_snp_file) as f:
        for line in f.readlines():
                if line.strip():
                        value=line.strip().split()
			snp[int(value[0])] = value[1]

called_snp_num, called_db_snp_num, called_match_db_snp_num = 0, 0, 0
called_indel_num, called_db_indel_num, called_match_db_indel_num = 0, 0, 0

called_snp_num = len(snp)
for key, value in snp.iteritems():
        if len(value) > 1:
                called_indel_num += 1
        if key-1 in db_snp:
                called_db_snp_num += 1
                if len(value) > 1:
                        called_db_indel_num += 1
                if value in db_snp[key-1]:
                        called_match_db_snp_num += 1
                        if len(value) > 1:
                                called_match_db_indel_num += 1
                                        
result_file_path = result_path + "/SRR352199-chr20-gatk-eval.txt"
result_file = open(result_file_path, "w")

result_file.write("SRR352199_1.fastq + SRR352199_2.fastq\n")
result_file.write("db_snp_num\t" + str(len(db_snp)) + "\n")
result_file.write("db_indel_num\t" + str(len(db_indel)) + "\n")
result_file.write("called_snp_num\t" + str(called_snp_num) + "\n")
result_file.write("called_indel_num\t" + str(called_indel_num) + "\n")
result_file.write("called_db_snp_num\t" + str(called_db_snp_num) + "\n")
result_file.write("called_db_indel_num\t" + str(called_db_indel_num) + "\n")
result_file.write("called_match_db_snp_num\t" + str(called_match_db_snp_num) + "\n")
result_file.write("called_match_db_indel_num\t" + str(called_match_db_indel_num) + "\n")

result_file.write("Called SNPs/db_snp_num\t" + str(float(called_db_snp_num)/len(db_snp)) + "\n")
result_file.write("Called Indels/db_indel_num\t" + str(float(called_db_indel_num)/len(db_indel)) + "\n")

result_file.close()
