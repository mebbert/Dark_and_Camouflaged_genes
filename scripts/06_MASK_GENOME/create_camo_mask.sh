#!/bin/bash

function maskGenome {
	expanded=$1
	i=$2
	REF=$3
	RESULT_DIR=$4
	PREFIX=$5

	genome=${REF}.fai
	ploidy=$((2*i))
	out="${RESULT_DIR}/${PREFIX}.ploidy_${ploidy}.fa"
	echo " awk "\$NF == $i" $expanded | bedtools complement -i - -g $genome | bedtools maskfasta -fi $REF -bed - -fo $out"
	awk "\$NF == $i" $expanded | \
		bedtools complement \
			-i - \
			-g $genome | \
		bedtools maskfasta \
			-fi $REF \
			-bed - \
			-fo $out

	## Prepare the new masked file for bwt
	echo " bwa index -a bwtsw b37-camo_mask.fa"
	bwa index -a bwtsw $out
	echo " samtools faidx b37-camo_mask.fa"
	samtools faidx $out
	echo "picard CreateSequenceDictionary REFERENCE=hg38-camo_mask.fa  OUTPUT=hg38-camo_mask.dict"
	picard CreateSequenceDictionary \
		REFERENCE=$out
		OUTPUT=${out%.*}.dict
}

camo_bed=$1
REF=$2
RESULT_DIR=$3
PREFIX=$4

mkdir -p $RESULT_DIR

echo " bedtools slop -b 50 -i $camo_bed -g $genome |  bedtools complement -i - -g $genome | bedtools maskfasta -fi $REF -bed - -fo b37-camo_mask.fa"
expanded=${camo_bed//.bed/.expanded_50.bed}
bedtools slop \
	-b 50 \
	-i $camo_bed \
	-g ${REF}.fai \
	> $expanded

for i in $(seq 2 1 5)
do
	maskGenome \
		$expanded \
		$i \
		$REF \
		$RESULT_DIR \
		$PREFIX &
done
wait

