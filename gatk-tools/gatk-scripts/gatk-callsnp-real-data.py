import os

ref_path = "/backup/SNP/GRCh37_human/genomes"
sam_path = "/data/nsvo/Human/SRS000204/results"
results_path = "/data/nsvo/Human/SRS000204/results"
script_path = "/data/nsvo/gatk-scripts"
tool_path = "/data/nsvo/genome-tools"

cmd = "time " + script_path + "/gatk-callsnp.sh " + ref_path + "/GRCh37_chr1.fasta " \
	+ sam_path + "/SRR352199-chr1_sorted_RG.bam " \
	+ ref_path + "/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes-diffcontigname.vcf " \
	+ results_path + "/SRR352199-chr1.vcf " + tool_path
os.system(cmd)
