echo "Calling SNPs using GATK"

echo 'Input FASTA file: ' $1
echo 'Input BAM files: ' $2
echo 'Input VCF files: ' $3
echo 'Output VCF files: ' $4
echo 'genome-tools path: ' $5

echo 'Start calling SNPs...'
java -jar $5/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO -nt 32 \
	-R $1 -T UnifiedGenotyper -I $2 --dbsnp $3 -o $4 \
	--output_mode EMIT_VARIANTS_ONLY -glm BOTH -rf BadCigar \
	-stand_call_conf 10.0 -stand_emit_conf 10.0 -dcov 200
echo 'Finish calling SNPs.'
