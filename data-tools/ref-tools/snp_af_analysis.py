import sys

snp_index_file = sys.argv[1]
epsilon = float(sys.argv[2])

f = open(snp_index_file)
count_epsilon, count_zero = 0, 0
for line in f:
	if line.strip():
		tmp = line.strip().split()
		snp_prof_len = (len(tmp) - 1)/2
		for i in range(snp_prof_len):
			if float(tmp[2*i+2]) <= epsilon:
				print tmp[2*i+2]
				count_epsilon += 1
				if float(tmp[2*i+2]) == 0:
					count_zero += 1
#print count_epsilon, count_zero