echo "Starting preprocessing data for GATK-UGT"

echo 'Input FASTA file (remove file extension .fasta): ' $1
echo 'Input SAM files (remove file extension .sam): ' $2
echo 'genome-tools path: ' $3

echo 'Sorting the BAM file...'
echo 'java -Xmx2g -Djava.io.tmpdir=/data/nsvo/tmp -jar $3/picard-tools-1.109/picard-tools-1.109/SortSam.jar \
	INPUT=$2.bam OUTPUT=$2_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/data/nsvo/tmp'
java -Xmx2g -Djava.io.tmpdir=/data/nsvo/tmp -jar $3/picard-tools-1.109/picard-tools-1.109/SortSam.jar \
	INPUT=$2.bam OUTPUT=$2_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/data/nsvo/tmp
echo 'Finish sorting the BAM file'

echo 'Creating BAM index file...'
echo '$3/samtools-0.1.19/samtools index $2_sorted.bam $2_sorted.bam.bai'
$3/samtools-0.1.19/samtools index $2_sorted.bam $2_sorted.bam.bai
echo 'Finish creating BAM index file'

echo 'Checking format requirement...'
echo "$3/samtools-0.1.19/samtools view -H $2_sorted.bam | grep='@RG'"
$3/samtools-0.1.19/samtools view -H $2_sorted.bam | grep='@RG'
echo 'Finish checking format requirement'

echo 'Adding read groups...'
echo "java -jar $3/picard-tools-1.109/picard-tools-1.109/AddOrReplaceReadGroups.jar \
	I=$2_sorted.bam O=$2_sorted_RG.bam SORT_ORDER=coordinate RGID=foo \
	RGLB=bar RGPL=illumina RGPU=run RGSM=DePristo CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT"
java -jar $3/picard-tools-1.109/picard-tools-1.109/AddOrReplaceReadGroups.jar \
	I=$2_sorted.bam O=$2_sorted_RG.bam SORT_ORDER=coordinate RGID=foo \
	RGLB=bar RGPL=illumina RGPU=run RGSM=DePristo CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT
echo 'Finish Adding read groups'

echo "Finish preprocesssing."
