import os
import sys

chr_id = sys.argv[1]
idx_dir = sys.argv[2]
cmd = "time go run /home/nsvo/workspace/goprojects/src/github.com/namsyvo/ISC/main/index.go \
	-g /data/nsvo/Human/refs/GRCh37_chr" + chr_id + "/GRCh37_chr" + chr_id + ".fasta \
	-s /data/nsvo/Human/refs/GRCh37_chr" + chr_id + "/ALL.chr" + chr_id + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf \
	-i /data/nsvo/test_data/GRCh37_chr" + chr_id + "/" + idx_dir + \
	" 2>/data/nsvo/test_data/GRCh37_chr" + chr_id + "/" + idx_dir + "/index.log"
os.system(cmd)
