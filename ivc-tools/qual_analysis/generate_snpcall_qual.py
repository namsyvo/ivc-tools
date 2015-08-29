import math

epsilon = 0.00001
Prior_prob = {}

#QualtoProb converts base qualities decoded by ASCII codes to probabilities
def QualtoProb(e):
	return math.pow(10, -(ord(e) - 33)/10.0)

#ProbtoQual converts probabilities to phred-scale quality scores
def ProbtoQual(p):
	return -10*math.log(1 - p, 10)

def CalcPostProb_ascii(snp):
	a = snp[0]
	q = snp[1]

	p = 0.0
	p_ab = {}
	p_a = 0.0

	for b, p_b in Prior_prob.iteritems():
		if a == b:
			p = 1.0 - math.pow(10, -(ord(q) - 33) / 10.0) #Phred-encoding factor (33)
		else:
			p = math.pow(10, -(ord(q) - 33) / 10.0) / 3
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	for b, p_b in Prior_prob.iteritems():
		Prior_prob[b] = p_b * (p_ab[b] / p_a)

def CalcPostProb(snp):
	a = snp[0]
	q = snp[1]

	p = 0.0
	p_ab = {}
	p_a = 0.0

	for b, p_b in Prior_prob.iteritems():
		if a == b:
			p = 1.0 - q
		else:
			p = q / 3
		p_ab[b] = p
		p_a += p_b * p_ab[b]
	for b, p_b in Prior_prob.iteritems():
		Prior_prob[b] = p_b * (p_ab[b] / p_a)

def InitProb():
	Prior_prob["A"] = 1 - 3*epsilon
	#Prior_prob["A"] = 0.5 - epsilon
	Prior_prob["C"] = epsilon
	#Prior_prob["C"] = 0.5 - epsilon
	Prior_prob["G"] = epsilon
	Prior_prob["T"] = epsilon

if __name__ == "__main__":

	for snp, prob in Prior_prob.iteritems():
		print "SNP: ", snp, "\t"
		print "Prob: ", prob, "\tQual: ", ProbtoQual(prob)

	'''
	snp_qual_ascii = [
		["C", "?"],
		#["G", "I"],
	]

	print "ASCII:"
	for snp in snp_qual_ascii:
		print "a: ", snp[0]
		print "q: ", snp[1]
		CalcPostProb_ascii(snp)
		for snp, prob in Prior_prob.iteritems():
			print "SNP: ", snp, "\t"
			print "Prob: ", prob, "\tQual: ", ProbtoQual(prob)

	snp_qual_prob = [
		#["A", 0.001],
		#["G", 0.000305],
		#["T", 0.00015],
		["G", 0.0002],
	]

	print "Prob:"
	for snp in snp_qual_prob:
		print "a: ", snp[0]
		print "q: ", snp[1]
		CalcPostProb(snp)
		for snp, prob in Prior_prob.iteritems():
			print "SNP: ", snp, "\t"
			print "Prob: ", prob, "\tQual: ", ProbtoQual(prob)
	'''

	snp_qual_prob_arr_A = [["A", 0.00005*i] for i in range(1, 51)]
	snp_qual_prob_arr_T = [["T", 0.00005*i] for i in range(1, 51)]
	f = open("prior_err_posterior.txt", "w")
	for i in range(len(snp_qual_prob_arr_A)):
		InitProb()
		f.write(str(snp_qual_prob_arr_T[i][1]) + "\t")
		#CalcPostProb(snp_qual_prob_arr_C[i])
		#CalcPostProb(snp_qual_prob_arr_C[i])
		#CalcPostProb(snp_qual_prob_arr_T[i])
		CalcPostProb(snp_qual_prob_arr_A[i])
		CalcPostProb(snp_qual_prob_arr_T[i])
		#f.write(str(Prior_prob[snp_qual_prob_arr_C[i][0]]) + "\t" + str(ProbtoQual(Prior_prob[snp_qual_prob_arr_C[i][0]])) + "\n")
		f.write(str(Prior_prob[snp_qual_prob_arr_T[i][0]]) + "\t" + str(ProbtoQual(Prior_prob[snp_qual_prob_arr_T[i][0]])) + "\n")
	f.close()
