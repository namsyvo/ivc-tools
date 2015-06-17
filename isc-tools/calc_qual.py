import math
import sys

p_q = sys.argv[1]
base_qual = sys.argv[2]

if p_q == "q2p":
	print -math.log(1 - math.pow(10, -(ord(base_qual) - 33)/10.0), 10)
elif p_q == "p2q":
	print -10 * math.log(1 - float(base_qual), 10)