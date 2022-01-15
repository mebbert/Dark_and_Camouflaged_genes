#!/bin/usr bash

# Enable bash debugging to log all commands
set -x

# Mask all but one of the identical regions for each given set of
# identical regions. i.e., if there are three identical regions, mask two
# of the three so all reads will align to a single region (to remove alignment
# ambiguity).
function maskGenome {
	EXPANDED_BED=$1
	UNMASKED_REF=$2
	PREFIX=$3

	unmasked_ref_index=${UNMASKED_REF}.fai
	masked_ref="${PREFIX}.fa"
	cat $EXPANDED_BED | \
		bedtools complement \
			-i - \
			-g $unmasked_ref_index | \
		bedtools maskfasta \
			-fi $UNMASKED_REF \
			-bed - \
			-fo $masked_ref

	## Prepare the new masked file for bwt
	bwa index -a bwtsw $masked_ref
	samtools faidx $masked_ref
	gatk CreateSequenceDictionary -R $masked_ref
}

ALIGN_TO_BED=$1
REF=$2
PREFIX=$3

# Create an expanded bed files to allow reads that overlap with the region to
# align easily (including on the ends of the region).
EXPANDED_BED=${ALIGN_TO_BED//.bed/.expanded_50.bed}
bedtools slop \
	-b 50 \
	-i $ALIGN_TO_BED \
	-g ${REF}.fai \
	> $EXPANDED_BED


maskGenome \
	$EXPANDED_BED \
	$REF \
	$PREFIX

