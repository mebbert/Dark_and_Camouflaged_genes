#!/bin/bash


cram=$1
ref=$2
ref_tag=$3
cram_ref=$4

#mkdir -p $out

echo "`date` Realign $1"
echo "`date` job running on `hostname`"

regex="SM:([A-Za-z0-9_\-]+)"
RG=$(samtools view -H $cram | grep '^@RG' | tail -1)
[[ $RG =~ $regex ]] 
sampleName=${BASH_REMATCH[1]}
RG=${RG//	/\\t}
echo "`date` Realigning sample: $sampleName" 
TMP_DIR="tmp/$sampleName"
mkdir -p $TMP_DIR

sorted_cram="${TMP_DIR}/${sampleName}.name_sorted.cram"

echo "`date` Sorting bam by name"
echo "`date` samtools sort -@ 16 -n -m 16G $cram > $sorted_cram"
if ! samtools sort -@ 16 -n -m 14G $cram > $sorted_cram; then
	echo "sort failed"
	exit 1
fi

fq1="${TMP_DIR}/${sampleName}_R1.fastq"
fq2="${TMP_DIR}/${sampleName}_R2.fastq"

#TODO - set the correct reference
export CRAM_REFERENCE="${cram_ref}"

echo "`date` Converting to Fastq:"
echo "`date` bedtools bamtofastq -i $sorted_cram -fq $fq1" # -fq2 $fq2"
if ! bedtools bamtofastq -i $sorted_cram -fq $fq1 -fq2 $fq2; then
	echo "bamtofastq failed"
	exit 1
fi

tmp_bam="${TMP_DIR}/${sampleName}.unsorted.bam"
echo "${ref}"

echo "`date` Aligning to GRCh38:"
echo "`date` bwa mem -M -R $RG -t 16 $ref $fq1 | samtools view -hb > $tmp_bam"
if ! bwa mem -M -R "$RG" -t 16 $ref $fq1 $fq2 | samtools view -hb > $tmp_bam; then
	echo "mem and veiw failed"
	exit 1
fi

final_bam="${sampleName}.${ref_tag}.bam"

echo "`date` Sorting final bam:"
echo "`date` samtools sort -@ 16 -m 16G $tmp_bam > $final_bam"
if ! samtools sort -@ 16 -m 14G $tmp_bam > $final_bam; then
	echo "sorting final bam failed"
	exit 1
fi

echo "`date` Indexing final bam:"
echo "`date` samtools index $final_bam"
if ! samtools index $final_bam; then
	echo "indexing final bam failed"
	exit 1
fi


echo "`date` Validating bam:"
if ! samtools quickcheck $final_bam; then
	echo "quickcheck failed"
	exit 1
fi


echo "`date` DONE!"
#rm -rfv $TMP_DIR
