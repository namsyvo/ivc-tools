import sys

dir_name = "/data/nsvo/test-data/GRCh37_chr1/reads/sim-reads/af_sid_mutant_dwgsim"
fn = sys.argv[1]
read_id = sys.argv[2]
inf = open(dir_name + "/" + fn)
outf = open(dir_name + "/alignment-analysis/" + fn + "." + read_id, "w")
while True:
    line = inf.readline()
    info = line.strip().split('_')
    if info[len(info)-1].split('/')[0] == read_id:
        outf.write(line)
        line = inf.readline()
        outf.write(line)
        line = inf.readline()
        outf.write(line)
        line = inf.readline()
        outf.write(line)
        break
outf.close()
