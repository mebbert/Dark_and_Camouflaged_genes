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
	if ! cat $EXPANDED_BED | \
		bedtools complement \
			-i - \
			-g $unmasked_ref_index | \
		bedtools maskfasta \
			-fi $UNMASKED_REF \
			-bed - \
			-fo $masked_ref; then
		echo "bedtools failed to mask the ref fasta"
		exit 1
	fi

	## Prepare the new masked file for bwt
	bwa index -a bwtsw $masked_ref
	samtools faidx $masked_ref
	gatk CreateSequenceDictionary -R $masked_ref
}

# This is called the 'align_to_bed' because it's used to mask the genome, leaving only a single
# 'camouflaged' region for each set of regions that are identical to each other.
ALIGN_TO_BED=$1
ALIGN_TO_REF=$2

# Prefix for output file names. Should be very specific to the reference genome version.
PREFIX=$3

# Create an expanded bed files to allow reads that overlap with the region to
# align easily (including on the ends of the region).
EXPANDED_BED=${ALIGN_TO_BED//.bed/.expanded_50.bed}
if ! bedtools slop \
	-b 50 \
	-i $ALIGN_TO_BED \
	-g ${ALIGN_TO_REF}.fai \
	> $EXPANDED_BED; then
	echo "bedtools failed to invert the bed file for masking"
	exit 1
fi


maskGenome \
	$EXPANDED_BED \
	$ALIGN_TO_REF \
	$PREFIX

