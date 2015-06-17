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
called_snp_file = result_path + "/SRR352199-chr1.vcf.txt"
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
                                        
result_file_path = result_path + "/SRR352199-chr1-gatk-eval.txt"
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

'''
if called_snp_num != 0 and snp_num != 0:
        result_file.write("SNP-Precision\t" + str(float(true_called_snp_num)/float(called_snp_num)) + "\n")
        result_file.write("SNP-Recall\t" + str(float(true_called_snp_num)/float(snp_num - called_snp_num + true_called_snp_num)) + "\n")
else:
        result_file.write("SNP-Precision\t\n")
        result_file.write("SNP-Recall\t\n")
        
if called_indel_num != 0 and indel_num != 0:
        result_file.write("INDEL-Precision\t" + str(float(true_called_indel_num)/float(called_indel_num)) + "\n")
        result_file.write("INDEL-Recall\t" + str(float(true_called_indel_num)/float(indel_num - called_indel_num + true_called_indel_num)) + "\n")
else:
        result_file.write("INDEL-Precision\t\n")
        result_file.write("INDEL-Recall\t\n")
        result_file.write("\n")
'''

result_file.close()

'''
log_file_path = result_path + "/SRR062641-log.txt"
log_file = open(log_file_path, "w")

stat_file_path = result_path + "/SRR062641-stat.txt"
stat_file = open(stat_file_path, "w")

with open(log_file) as f:
        lines = f.readlines()

        stat_file.write("SRR062641_1.fastq\n")
        stat_file.write(lines[2])
        stat_file.write(lines[3])
        stat_file.write(lines[4])
        stat_file.write(lines[6])
        stat_file.write(lines[7])
        
stat_file.close()
'''
