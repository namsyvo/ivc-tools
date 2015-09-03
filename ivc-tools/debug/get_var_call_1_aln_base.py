"""
Get var_calls with 1_base_aligned
"""
import sys
if len(sys.argv) != 3:
    print "Usage: python get_aln_info_diff_var.py var_type result_file"
    exit(0)

var_type = sys.argv[1]
result_file = sys.argv[2]

var_pos = {}
file_out = open(result_file + ".1base_aln", "w")
file_in = open(result_file)
line = file_in.readline()
tmp = line.strip().split()
if var_type == "tp":
    file_out.write("\t".join(tmp[0:19]) + "\ttpfp\n")
elif var_type == "fp":
    file_out.write("\t".join(tmp[0:1] + tmp[3:21]) + "\ttpfp\n")

for line in file_in.readlines():
    tmp = line.strip().split()
    #problem of input data: colum base_qual is duplicate
    if var_type == "tp" and int(tmp[7]) == 1:
        file_out.write("\t".join(tmp[0:9] + tmp[10:20]) + "\ttpfp\n")
    elif var_type == "fp" and int(tmp[9]) == 1:
        file_out.write("\t".join(tmp[0:1] + tmp[3:11] + tmp[12:22]) + "\ttpfp\n")
file_out.close()