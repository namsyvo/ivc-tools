echo "Starting preprocessing data for GATK-UGT"

echo 'Input FASTA file (remove file extension .fasta): ' $1
echo 'Input SAM files (remove file extension .sam): ' $2
echo 'genome-tools path: ' $3

echo 'Converting SAM file to BAM file...'
$3/samtools-0.1.19/samtools view -bS -t $1.fasta $2.sam > $2.bam
echo 'Finish converting SAM file to BAM file'

echo 'Creating reference dictionary...'
java -jar $3/picard-tools-1.109/picard-tools-1.109/CreateSequenceDictionary.jar \
	R=$1.fasta O=$1.dict
echo 'Finish creating reference dictionary'

echo 'Creating fasta index file...'
$3/samtools-0.1.19/samtools faidx $1.fasta
echo 'Finish creating fasta index file'

echo 'Checking format requirement...'
$3/samtools-0.1.19/samtools view -H $2.bam
echo 'Finish checking format requirement'

echo 'Sorting the BAM file...'
echo 'java -jar $3/picard-tools-1.109/picard-tools-1.109/SortSam.jar \
	INPUT=$2.bam OUTPUT=$2_sorted.bam SORT_ORDER=coordinate'
java -jar $3/picard-tools-1.109/picard-tools-1.109/SortSam.jar \
	INPUT=$2.bam OUTPUT=$2_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
echo 'Finish sorting the BAM file'

#echo 'Resorting the BAM file...'
#java -jar $3/picard-tools-1.109/ReorderSam.jar \
#	INPUT=$2.bam OUTPUT=$2_resorted.bam
#echo 'Finish resorting the BAM file'

echo 'Creating BAM index file...'
$3/samtools-0.1.19/samtools index $2_sorted.bam $2_sorted.bam.bai
echo 'Finish creating BAM index file'

echo 'Checking format requirement...'
$3/samtools-0.1.19/samtools view -H $2_sorted.bam | grep='@RG'
echo 'Finish checking format requirement'

echo 'Adding read groups...'
java -jar $3/picard-tools-1.109/picard-tools-1.109/AddOrReplaceReadGroups.jar \
	I=$2_sorted.bam O=$2_sorted_RG.bam SORT_ORDER=coordinate RGID=foo \
	RGLB=bar RGPL=illumina RGPU=run RGSM=DePristo CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT
echo 'Finish Adding read groups'

echo "Finish preprocesssing."
