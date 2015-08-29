import math
import sys

p_q = sys.argv[1]
base_qual = sys.argv[2]

if p_q == "q2p":
	print ord(base_qual), math.pow(10, -(ord(base_qual) - 33)/10.0)
elif p_q == "p2q":
	q_num = -10 * math.log(float(base_qual), 10)
	print str(unichr(int(32+q_num))), q_num