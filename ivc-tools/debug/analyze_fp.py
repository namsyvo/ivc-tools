import sys
import math

f = open(sys.argv[1])
data = f.readlines()
mis = 0
total = len(data)
for line in data:
	info = line.split()
	if info[6] != 'None' and info[9] != 'None':
		if math.fabs(int(info[6])) == math.fabs(int(info[9])):
			mis += 1
print mis
print total