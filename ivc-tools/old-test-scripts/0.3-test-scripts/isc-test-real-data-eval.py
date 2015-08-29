import sys

confi = float(sys.argv[1])
snp_file = sys.argv[2]
result_path = sys.argv[3]
time = sys.argv[4]

db_snp = {}
db_indel = {}

with open(snp_file) as f:
        for line in f.readlines():
                if line.strip():
                        value=line.strip().split()
                        db_snp[int(value[0])] = value[1:]
                        if len(value[1]) > 1:
                                db_indel[int(value[0])] = value[1:]

snp = {}
called_snp_file = result_path + "/SRR062641-full-32-called-snps.txt"
with open(called_snp_file) as f:
        for line in f.readlines():
                if line.strip():
                        value=line.strip().split()
                        if len(value) == 3:
                                if float(value[2]) >= confi:
                                        snp[int(value[0])] = value[1]

called_snp, called_known_allele_snp, called_known_snp = 0, 0, 0
called_indel, called_known_allele_indel, called_known_indel = 0, 0, 0

called_snp = len(snp)
for key, value in snp.iteritems():
        if len(value) > 1:
                called_indel += 1
        if key in db_snp:
                called_known_snp += 1
                if len(value) > 1:
                        called_known_indel += 1
                if value in db_snp[key]:
                        called_known_allele_snp += 1
                        if len(value) > 1:
                                called_known_allele_indel += 1

result_file_path = result_path + "/snp-num-time-mem-" + str(confi) + ".txt"
result_file = open(result_file_path, "w")

result_file.write("Alg\trun\tread\tsize\tproc\ttimeI\tmemI\ttimeA\tmemA\tsnp\tindel\tcalled_snp\tcalled_known_snp\tcalled_known_allele_snp\tcalled_indel\tcalled_known_indel\tcalled_known_allele_indel\tknown_allele_snp_ratio\tknown_allele_indel_ratio\n")
result_file.write("isc-0.3.6-snpqual-inalignment_5168d045aaf4d7db5d872f798be73d1ceb026a00+fmisnp-gr_797620715c769a731d449179f62f92af2a185b1d\t" + time + "\tSRR062641_1.fastq-SRR062641_2.fastq\tfull\t32\t")

mem_time_file = result_path + "/SRR062641-full-32-called-snps.txt.log"
f = open(mem_time_file)
for line in f:
    tokens = line.strip().split("\t")
    if "time for initializing SNP caller" in tokens[0]:
        result_file.write(tokens[1] + "\t")
    if "memstats after initializing SNP caller" in tokens[0]:
        result_file.write(tokens[3] + "\t")
    if "time for calling SNPs" in tokens[0]:
        result_file.write(tokens[1] + "\t")
    if "memstats after calling SNPs" in tokens[0]:
        result_file.write(tokens[3] + "\t")
f.close()

result_file.write(str(len(db_snp)) + "\t" + str(len(db_indel)) + "\t")
result_file.write(str(called_snp) + "\t" + str(called_known_snp) + "\t" + str(called_known_allele_snp) + "\t")
result_file.write(str(called_indel) + "\t" + str(called_known_indel) + "\t" + str(called_known_allele_indel) + "\t")
result_file.write(str(float(called_known_allele_snp)/len(db_snp)) + "\t" + str(float(called_known_allele_indel)/len(db_indel)) + "\n")

result_file.close()
