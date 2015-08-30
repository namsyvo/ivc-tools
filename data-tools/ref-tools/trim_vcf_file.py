import sys
import os

vcf_file = open(sys.argv[1])
dir_name, file_name = os.path.split(sys.argv[1])
trim_vcf_file = open(os.path.join(dir_name, "TRIMMED." + file_name), "w")

line = vcf_file.readline()
while True:
    if line[:2] != "##":
        break
    trim_vcf_file.write(line)
    line = vcf_file.readline()

while line:
    trim_vcf_file.write("\t".join(line.split()[:9]) + "\n")
    line = vcf_file.readline()

trim_vcf_file.close()
