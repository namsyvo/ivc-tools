#Preprocess data for GATK-UGT:
python bwa-align-running.py.py
python gatk-preprocess-running.py
go run vcf_replace_chr.go vcf_file > vcf_rep_contig

#Run GATK and evalate results:
python gatk-callsnp-running.py
python extract-alt-vcf-running.py
python gatk-eval.py
