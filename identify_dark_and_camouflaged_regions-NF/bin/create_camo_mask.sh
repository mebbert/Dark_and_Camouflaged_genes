#!/bin/usr bash

# Enable bash debugging to log all commands
set -x

# Mask all but one of the identical regions for each given set of
# identical regions. i.e., if there are three identical regions, mask two
# of the three so all reads will align to a single region (to remove alignment
# ambiguity).
function maskGenome {
	expanded=$1
	i=$2
	REF=$3
	PREFIX=$4

	genome=${REF}.fai
	ploidy=$((2*i))
	out="${PREFIX}.ploidy_${ploidy}.fa"
	awk "\$NF == $i" $expanded | \
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

camo_bed=$1
REF=$2
PREFIX=$3

# Create an expanded bed files to allow reads that overlap with the region to
# align easily (including on the ends of the region).
expanded=${camo_bed//.bed/.expanded_50.bed}
bedtools slop \
	-b 50 \
	-i $camo_bed \
	-g ${REF}.fai \
	> $expanded


# For set of regions with identical matches from 2:5, mask the genome. We deal
# with each ploidy individually because it affects how we call variants. Ploidy
# goes up by the number of regions that are identical. e.g., if there are three
# regions that are identical, ploidy is 2*3 = 6. We stop at five identical regions
# because we do not have much confidence in calls beyond that; it's an area of
# research.
for i in $(seq 2 1 5)
do
	maskGenome \
		$expanded \
		$i \
		$REF \
		$PREFIX &
done

# Wait for all jobs launched in background to complete. Then return.
wait

