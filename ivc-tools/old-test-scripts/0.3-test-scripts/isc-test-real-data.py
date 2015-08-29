import os
import sys
from datetime import datetime

config_fn = sys.argv[1]
f = open(config_fn)
prog_version = f.readline().strip()
prog_path = f.readline().strip()
ref_path = f.readline().strip()
genome_fn = f.readline().strip()
snp_fn = f.readline().strip()
read_path = f.readline().strip()
read_fn = f.readline().strip()
data_path = f.readline().strip()
index_dn = f.readline().strip()
f.close()

confi = float(sys.argv[2])
cpu_num = int(sys.argv[3])
cpu_use = sys.argv[4]
sizes = sys.argv[5:]
proc_nums = []
if cpu_use == "only":
        proc_nums = [str(cpu_num)]
elif cpu_use == "all":
        for i in range(1, 1 + cpu_num/2):
                proc_nums.append(str(2*i))

genome_file = ref_path + "/" + genome_fn
snp_file = ref_path + "/" + snp_fn

index_dir = data_path + "/" + index_dn
time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "-")
result_path = data_path + "/results/real-reads/" + prog_version + "-" + time_stamp
if not os.path.exists(result_path):
    os.makedirs(result_path)

for size in sizes:
        for proc_num in proc_nums:
                query_file_1 = read_path + "/" + read_fn + "_1-" + size + ".fastq"
                query_file_2 = read_path + "/" + read_fn + "_2-" + size + ".fastq"

                called_snp_file = result_path + "/" + read_fn + "-" + size + "-" + proc_num + "-called-snps.txt"
                mem_time_file = result_path + "/" + read_fn + "-" + size + "-" + proc_num + "-called-snps.txt.log"
                cmd = "(go run " + prog_path + \
                    " -g " + genome_file + " -s " + snp_file + " -i " + index_dir + \
                    " -1 " + query_file_1 + " -2 " + query_file_2 + " -o " + called_snp_file + \
                    " -w " + proc_num + " -t " + proc_num + ") 2>" + mem_time_file            
                print cmd
                os.system(cmd)

db_snp = {}
db_indel = {}

with open(index_dir + "/" + snp_fn + ".idx") as f:
        for line in f.readlines():
                if line.strip():
                        value = line.strip().split()
                        db_snp[int(value[0])] = value[1:]
                        if len(value[1]) > 1:
                                db_indel[int(value[0])] = value[1:]
result_file_path = result_path + "/snp-num-time-mem-" + str(confi) + ".txt"
result_file = open(result_file_path, "w")
result_file.write("Alg\trun\tread\tsize\tproc\ttimeI\tmemI\ttimeA\tmemA\tsnp\tindel\tcalled_snp\tcalled_known_snp\tcalled_known_allele_snp\tcalled_indel\tcalled_known_indel\tcalled_known_allele_indel\tknown_allele_snp_ratio\tknown_allele_indel_ratio\n")

for size in sizes:
        for proc_num in proc_nums:

                result_file.write(prog_version + "\t" + time_stamp + "\t" + read_fn + "\t" + size + "\t" + proc_num + "\t")

                mem_time_file = result_path + "/" + read_fn + "-" + size + "-" + proc_num + "-called-snps.txt.log"
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

                called_snp_file = result_path + "/" + read_fn + "-" + size + "-" + proc_num + "-called-snps.txt"
                snp = {}
                with open(called_snp_file) as f:
                        for line in f.readlines():
                                if line.strip():
                                        value=line.strip().split()
                                        if len(value) == 3:
                                                if float(value[2]) >= confi:
                                                        snp[int(value[0])] = value[1]

                snp_num, indel_num = len(db_snp), len(db_indel)
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

                result_file.write(str(snp_num) + "\t" + str(indel_num) + "\t")
                result_file.write(str(called_snp) + "\t" + str(called_known_snp) + "\t" + str(called_known_allele_snp) + "\t")
                result_file.write(str(called_indel) + "\t" + str(called_known_indel) + "\t" + str(called_known_allele_indel) + "\t")
                result_file.write(str(float(called_known_allele_snp)/len(db_snp)) + "\t" + str(float(called_known_allele_indel)/len(db_indel)) + "\n")

                result_file.close()
                print "Finish evaluating results for ", size, " ", proc_num