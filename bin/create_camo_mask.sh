#!/bin/usr bash

function maskGenome {
	expanded=$1
	i=$2
	REF=$3
	PREFIX=$4

	genome=${REF}.fai
	ploidy=$((2*i))
	out="${PREFIX}.ploidy_${ploidy}.fa"
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
	java -jar /usr/bin/picard/build/libs/picard.jar CreateSequenceDictionary \
		REFERENCE=$out
		OUTPUT=${out%.*}.dict
}

camo_bed=$1
REF=$2
PREFIX=$3

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
		$PREFIX &
done
wait

