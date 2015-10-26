import sys
vcf_file = open(sys.argv[1])
for line in vcf_file:
    if line[0] != '#':
        tokens = line.split()
        if len(tokens[3]) > 1 and len(tokens[4]) > 1:
            print line
