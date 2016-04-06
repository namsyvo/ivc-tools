
'''
#ref_path = "/home/SNP/NC_016088.1_soybean/genomes"
ref_path = "/home/SNP/GRCh37_human/genomes/"
ref = []
#f=open(ref_path + "/chr1.fasta")
#f=open(ref_path + "/GRCh37_chr20.fasta")
#f=open(ref_path + "/CM000682-chr20.fasta")
f=open(ref_path + "/human_g1k_v37.fasta")
f.readline()
for line in f.readlines():
	if line.strip():
		ref.extend(line.strip())
print len(ref)
print ref[0], ref[len(ref)-1]
'''
#ref_path = "/home/SNP/NC_016088.1_soybean/chr1_mutate1/index"
ref_path = "/home/SNP/GRCh37_human/index/"
ref = []
f=open(ref_path + "/genomestar.txt")
#f=open(ref_path + "/GRCh37_chr20.fasta")
for line in f.readlines():
	if line.strip():
		ref.extend(line.strip())
print len(ref)
print ref[0], ref[len(ref)-1]

#ref_path = "/home/SNP/NC_016088.1_soybean/genomes"
ref_path = "/home/SNP/GRCh37_human/genomes/"
#f=open(ref_path + "/vcf_chr_1.vcf")
#f=open(ref_path + "/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf")
#f=open(ref_path + "/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf")
with open(ref_path + "/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf") as f:
    for line in f:
		if line[0] != '#':
			#print line[0]
			info = line.strip().split()
			if ref[int(info[1])-1] != '*':
			#if ref[int(info[1])-1] != info[3][0]:
				print "*", info[1]
			'''
			if ref[int(info[1])-1] == 'N':
				print 'N', info[1]
			if info[2][0] != 'r' and info[2][0] != '.':
				print info[2]
			'''
