#!/bin/bash

## Here we take a reference genome and the Camo align_to file created
# In step 6 and mask the genome so only the align to bed regions are unmasked
# We only mask camo regions with repeat number ≤ 5, and have a seperate masked
# Genome for each ploidy we test

##################################################################
## TODO: Fill in correct paths                                  ##
##################################################################
ILL100_ALIGN_TO_BED="../../results/illuminaRL100/b37/illuminaRL100.b37.camo.align_to.sorted.bed" #Align_to camo bed created in step 5
ILL100_REALIGN_BED="../../results/illuminaRL100/b37/illuminaRL100.b37.camo.realign.sorted.bed" #Realign camo bed created in step 5
B38_ANNOTATIONS="../../results/annotation/b38/Homo_sapiens.GRCh38.93.annotation.bed" # b38 annotations created in step 4 (used to convert b37 beds to b38 coordinates)
B38_RESULT_DIR="../../results/illuminaRL100/b38" # where to write bed files with b38 coordinates


REF="Homo_saiens.GRCh37.dna.primary_assembly.sorted.fa" #b37 reference
MASK_RESULT_DIR="../../results/b37_camo_mask" #where to write masked genomes
MASK_REF_PREFIX="b37-camo_mask" #Prefix given to each of the masked genomes: (will be suffixed with .ploidy_X.fa depending on ploidy)
##################################################################

BASE=$(basename $ILL100_ALIGN_TO_BED)
B38_ALIGN_TO="${B38_RESULT_DIR}/${BASE//b37/b37}"
B38_REALIGN="${B38_ALIGN_TO//align_to/realign}"

python convert_to_hg38.py \
	$ILL100_ALIGN_TO_BED \
	$B38_ANNOTATIONS \
	> $B38_ALIGN_TO

python convert_to_hg38.py \
	$ILL100_REALIGN_BED \
	$B38_ANNOTATIONS \
	> $B38_REALIGN

bash create_camo_mask.sh \
	$ILL100_ALIGN_TO_BED \
	$REF \
	$MASK_RESULT_DIR \
	$MASK_REF_PREFIX
