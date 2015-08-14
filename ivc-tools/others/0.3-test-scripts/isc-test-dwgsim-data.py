import os
import sys
from datetime import datetime

config_file = sys.argv[1]
f=open(config_file)
prog_version = f.readline().strip()
prog_path = f.readline().strip()
data_path = f.readline().strip()
genome_fn = f.readline().strip()
snp_fn = f.readline().strip()
read_fn = f.readline().strip()
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

time_stamp = str(datetime.now())
time_stamp = time_stamp.replace(" ", "-")

ref_path = data_path + "/refs"
read_path = data_path + "/reads/sim-reads/dwgsim-data"
result_path = data_path + "/results/sim-reads/dwgsim-data/" + prog_version + "-" + time_stamp

if not os.path.exists(result_path):
    os.makedirs(result_path)

genome_file = ref_path + "/" + genome_fn
snp_file = ref_path + "/" + snp_fn
index_dir = data_path + "/" + index_dn

for size in sizes:
    for proc_num in proc_nums:
	query_file_1 = read_path + "/" + read_fn + "." + size + ".bwa.read1.fastq"
	query_file_2 = read_path + "/" + read_fn + "." + size + ".bwa.read2.fastq"

	called_snp_file = result_path + "/" + read_fn + "-" + size + "-" + proc_num + "-called-snps.txt"
	mem_time_file = result_path + "/" + read_fn + "-" + size + "-" + proc_num + "-called-snps.txt.log"

	cmd = "(go run " + prog_path + \
            " -g " + genome_file + " -s " + snp_file + " -i " + index_dir + \
            " -1 " + query_file_1 + " -2 " + query_file_2 + " -o " + called_snp_file + \
            " -w " + proc_num + " -t " + proc_num + ") 2>" + mem_time_file
	print cmd
	os.system(cmd)

true_snp = {}
true_indel = {}
#eval_file = data_path + "/refs/eval.txt"
eval_file = ref_path + "/GRCh37_chr20-mutate1/known_alleles.txt"
with open(eval_file) as f:
	for line in f.readlines():
		if line.strip():
			value = line.strip().split()
			true_snp[int(value[0])] = value[1]
			if len(value[1]) > 1:
				true_indel[int(value[0])] = value[1]

result_file_path = result_path + "/prec-rec-time-mem-" + str(confi) + ".txt"
result_file = open(result_file_path, "w")
result_file.write("Alg\trun\tread\tsize\tproc\ttimeI\tmemI\ttimeA\tmemA\tsnp\tindel\tcalled_snp\ttrue_called_snp\tcalled_indel\ttrue_called_indel\tsnp_prec\tsnp_rec\tindel_prec\tindel_rec\n")

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
                                    if len(value[2]) > 1 and float(value[2]) >= confi:
					snp[int(value[0])] = value[1]
	snp_num, indel_num = len(true_snp), len(true_indel)
	called_snp_num, true_called_snp_num = 0, 0
	called_indel_num, true_called_indel_num = 0, 0
	for key, value in snp.iteritems():
		if key in true_snp:
			called_snp_num += 1
			if len(value) > 1:
				called_indel_num += 1
			if value == true_snp[key]:
				true_called_snp_num += 1
				if len(value) > 1:
					true_called_indel_num += 1

	result_file.write(str(snp_num) + "\t" + str(indel_num) + "\t")
	result_file.write(str(called_snp_num) + "\t" + str(true_called_snp_num) + "\t")
	result_file.write(str(called_indel_num) + "\t" + str(true_called_indel_num) + "\t")
	if called_snp_num != 0 and snp_num != 0:
		result_file.write(str(float(true_called_snp_num)/float(called_snp_num)) + "\t")
		result_file.write(str(float(true_called_snp_num)/float(snp_num - called_snp_num + true_called_snp_num)) + "\t")
	else:
		result_file.write("\t\t")
	if called_indel_num != 0 and indel_num != 0:
		result_file.write(str(float(true_called_indel_num)/float(called_indel_num)) + "\t")
		result_file.write(str(float(true_called_indel_num)/float(indel_num - called_indel_num + true_called_indel_num)) + "\t")
	else:
		result_file.write("\t\t")
	result_file.write("\n")
result_file.close()
