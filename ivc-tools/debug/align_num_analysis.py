fn = "/data/nsvo/test-data/GRCh37_chr1/results/sim-reads/af_sid_mutant_dwgsim/ivc_0.70/IVC-0.5.4-aff_gap_aln-full_var-2015-06-16-22:29:23.571231/align-num-log.txt"
f = open(fn)
for j in range(0, 5):
        f.readline()
a = []
na = 0
for line in f:
        if "Outputing" in line:
                break
        if line.strip():
                a.append(int(line))
        else:
                na += 1
print len(a)
print na

#count = [[x, a.count(x)] for x in set(a)]
#print count

from collections import Counter
for key, value in Counter(a).items():
        print key, value
