iter_num_03, iter_num_031 = [], []
with open("real-test-500k-0.3.3-iter-num.txt") as f:
	for line in f.readlines():
		tmp = line.strip()
		iter_num_03.append(int(tmp))

print len(iter_num_03), sum(iter_num_03)
'''
with open("real-test-500k-0.3.1-iter-num.txt") as f:
        for line in f.readlines():
                tmp = line.strip()
                iter_num_031.append(int(tmp))

print len(iter_num_031), sum(iter_num_031)
'''
'''
iter_num_03, iter_num_031 = [], []
read_03, read_031 = [], []

with open("real-test-500k-0.3.txt") as f:
	for line in f.readlines():
		tmp = line.strip().split("\t")
		if int(tmp[0]) != 1:
			iter_num_03.append(int(tmp[0]))
			read_03.append(tmp[2])
with open("real-test-500k-0.3.1.txt") as f:
	for line in f.readlines():
		tmp = line.strip().split("\t")
		if int(tmp[0]) != 1:
			iter_num_031.append(int(tmp[0]))
			read_031.append(tmp[2])

print len(read_031), len(read_03)
diff = 0
for read in read_031:
	for i, r in enumerate(read_03):
	    if r == read and iter_num_031[i] == 6 and iter_num_03[i] < 6:
			print iter_num_03[i], iter_num_031[i], read

for read in read_031:
	for i, r in enumerate(read_03):
	    if r == read and iter_num_031[i] == 5 and iter_num_03[i] < 5:
			print iter_num_03[i], iter_num_031[i], read

for read in read_031:
	for i, r in enumerate(read_03):
	    if r == read and iter_num_031[i] == 4 and iter_num_03[i] < 4:
			print iter_num_03[i], iter_num_031[i], read
'''
