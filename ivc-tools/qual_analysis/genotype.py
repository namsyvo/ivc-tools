import math

def print_prof(prof):
	for k,v in prof.items():
		print("\t%s: %.3f" % (k, v))

def log2(x):
	return math.log(x) / math.log(2.0)

def H(prof):
	return sum( -prof[x] * log2(prof[x]) for x in list("ACGT") if prof[x]>0)

def prob_base_given_gen(a, q, g):
	if g.count(a) == 2:
		return 1-q
	elif g.count(a) == 1:
		return (3.0-2.0*q)/6.0
	else:
		return q/3

# new_prof[b] = prob that genotype is b, given aligned base is b
def P(b, q, prof):
	new_prof = {}
	s = sum(prof[g] * prob_base_given_gen(b, q, g) for g in Genotype)
	for g in Genotype:
		new_prof[g] = prof[g] * prob_base_given_gen(b, q, g) / s
	return new_prof

BASE_ACC = 0.999
Genotype = ["AA","AC","CC"]
prof = dict(AA=0.9985, AC=0.001, CC=0.0005)
#Genotype = ["AA","AC","AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"]
#prof = dict(AA=0.991, AC=0.001, AG=0.001, AT=0.001, CC=0.001, CG=0.001, CT=0.001, GG=0.001, GT=0.001, TT=0.001)

AlignedBases = ""
for i in range(63):
	AlignedBases += 'A'
for i in range(7):
	AlignedBases += 'C'
AlnBaseQual = []
for i in range(63):
	AlnBaseQual.append(35)
for i in range(7):
	AlnBaseQual.append(35)

print(prof)
for i in range(len(AlignedBases)):
	b = AlignedBases[i]
	q = math.pow(10, -AlnBaseQual[i]/10.0)
	prof = P(b, q, prof)
	print(b)
	for k,v in prof.items():
		print("\t%s: %.20f" % (k, v))
