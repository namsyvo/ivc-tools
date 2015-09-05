import sys
import os

vcf_file = open(sys.argv[1])
dir_name, file_name = os.path.split(sys.argv[1])
trim_vcf_file = open(os.path.join(dir_name, "TRIMMED." + file_name), "w")

for line in vcf_file:
    if line[:2] == "##":
        trim_vcf_file.write(line)
    else:
        trim_vcf_file.write("\t".join(line.split()[:9]) + "\n")
trim_vcf_file.close()
