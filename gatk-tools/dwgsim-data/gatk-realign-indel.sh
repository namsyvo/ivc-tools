echo "Calling SNPs using GATK"

echo 'Input FASTA file: ' $1
echo 'Input BAM files: ' $2
echo 'Input VCF files: ' $3
echo 'Output VCF files: ' $4
echo 'genome-tools path: ' $5

echo 'Start realigning indels...'
java -jar $5/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO -nt 32 \
	-T RealignerTargetCreator -R $1 -I $2 -o "$4.forIndelRealigner.intervals"

java -jar $5/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO \
	-T IndelRealigner -R $1 -I $2 -o $4 -targetIntervals "$4.forIndelRealigner.intervals"
echo 'Finish realigning indels.'
