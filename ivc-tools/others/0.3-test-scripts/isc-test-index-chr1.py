import os

cmd = "time go run /home/nsvo/workspace/goprojects/src/github.com/namsyvo/ISC/main/index.go \
	-g /data/nsvo/Human/refs/GRCh37_chr1/GRCh37_chr1.fasta \
	-s /data/nsvo/Human/refs/GRCh37_chr1/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf \
	-i /data/nsvo/test_data/GRCh37_chr1/index-fmisnp-gr-sep5 \
	2>/data/nsvo/test_data/GRCh37_chr1/index-fmisnp-gr-sep5/index.log"
os.system(cmd)
