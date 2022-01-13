#!/bin/usr bash

# Enable bash debugging to log all commands
set -x

# Mask all but one of the identical regions for each given set of
# identical regions. i.e., if there are three identical regions, mask two
# of the three so all reads will align to a single region (to remove alignment
# ambiguity).
function maskGenome {
	expanded_bed=$1
	REF=$2
	PREFIX=$3

	genome=${REF}.fai
	out="${PREFIX}.fa"
	awk "\$NF == $i" $expanded_bed | \
		bedtools complement \
			-i - \
			-g $genome | \
		bedtools maskfasta \
			-fi $REF \
			-bed - \
			-fo $out

	## Prepare the new masked file for bwt
	bwa index -a bwtsw $out
	samtools faidx $out
	java -jar /usr/bin/picard/build/libs/picard.jar CreateSequenceDictionary \
		REFERENCE=$out
		OUTPUT=${out%.*}.dict
}

align_to_bed=$1
REF=$2
PREFIX=$3

# Create an expanded bed files to allow reads that overlap with the region to
# align easily (including on the ends of the region).
expanded_bed=${align_to_bed//.bed/.expanded_50.bed}
bedtools slop \
	-b 50 \
	-i $align_to_bed \
	-g ${REF}.fai \
	> $expanded_bed


maskGenome \
	$expanded_bed \
	$REF \
	$PREFIX

