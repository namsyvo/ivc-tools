import sys
import math

f_in = open(sys.argv[1])
f_out = open(sys.argv[1] + ".qual", "w")
for line in f_in:
	tmp = line.split()
	tmp[5] = str(math.pow(10, -(ord(tmp[5][0]) - 33)/10.0))
	new_line = "\t".join(tmp[:14])
	f_out.write(new_line + "\n")
f_in.close()
f_out.close()