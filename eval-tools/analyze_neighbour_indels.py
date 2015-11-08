import argparse
parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str)
args = parser.parse_args()

stat = {}

with open(args.filename) as f:
   header = f.readline()
   for line in f:
      w, kw, t = line.strip().split("\t")
      w = int(w)
      kw = float(kw)
      if w not in stat:
         stat[w] = dict()
      if kw not in stat[w]:
         stat[w][kw] = dict(tp_snp=1, tp_indel=1, fn_snp=1, fn_indel=1)
      stat[w][kw][t] += 1

for w in stat:
   print(w)
   print("Variant\tKw\tTP\tFN\tTPFN")
   for k,v in sorted(stat[w].items()):
      print "SNP\t", round(k, 2), "\t", v['tp_snp'], "\t", v['fn_snp'], "\t", round(float(v['tp_snp'])/float(v['fn_snp']), 2)
   for k,v in sorted(stat[w].items()):
      print "INDEL\t", round(k, 2), "\t", v['tp_indel'], "\t", v['fn_indel'], "\t", round(float(v['tp_indel'])/float(v['fn_indel']), 2)
