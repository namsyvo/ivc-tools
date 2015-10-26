echo "Calling SNPs using GATK"

echo 'Input FASTA file: ' $1
echo 'Input BAM files: ' $2
echo 'Input VCF files: ' $3
echo 'Output VCF files: ' $4
echo 'genome-tools path: ' $5

echo 'Start calling variants...'
java -jar $5/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO -nct 32 \
	-T HaplotypeCaller -R $1 -I $2 -o $4 \
	--genotyping_mode DISCOVERY -rf BadCigar \
	-stand_call_conf 10.0 -stand_emit_conf 10.0 -dcov 200
echo 'Finish calling SNPs.'
