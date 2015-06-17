import os

result_path = "/data/nsvo/Human/SRS000204/results"
script_path = "/home/nsvo/workspace/isc-scripts/gatktools"

cmd = "go run " + script_path + "/extract-alt-vcf.go " + result_path + "/SRR352199-chr1.vcf > " \
	+ result_path + "/SRR352199-chr1.vcf.txt"
os.system(cmd)
