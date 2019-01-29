#!/bin/bash

# INTERSECT VCF
# intersect VCF with the camo region bed to find just variants that are within these homologous regions
# Should remove a lot of false positive variants in regions where reads were being forced to align

#Also since this is exome data restrict variants to only be found within CDS exons!! (We expect less
#false positives to be found in CDS exons

INPUT_VCF=$1
CAMO_BED=$2
ANNO_BED=$3
RESULT_DIR=$4
CDS_BED="/tmp/CDS_annos.bed"

awk '$4 == "CDS"' $ANNO_BED > $CDS_BED
bedtools intersect \
	-header \
	-a $INPUT_VCF \
	-b $CAMO_BED | \
	bedtools intersect \
	-header \
	-a - \
	-b $CDS_BED \
