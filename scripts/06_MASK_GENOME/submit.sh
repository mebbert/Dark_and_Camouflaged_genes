#!/bin/bash

## 06_MASK_GENOME

## In step 6, we take a reference genome and the Camo align_to file created
# In step 5 and mask the genome so only the align_to bed regions are unmasked
# We only mask camo regions with repeat number ≤ 5, and have a seperate masked
# Genome for each ploidy we test

##################################################################
## TODO: Fill in correct paths                                  ##
##################################################################
ILL100_ALIGN_TO_BED="../results/hg38_no_alt/illuminaRL100/illuminaRL100.hg38_no_alt.camo.align_to.sorted.bed" #Align_to camo bed created in step 5

REF="/research/labs/moleneurosci/petrucelli/m199515/NGS_resources/references/refdata-GRCh38-2.1.0/fasta/genome.fa" #b37 reference
MASK_RESULT_DIR="../results/hg38_camo_mask" #where to write masked genomes
MASK_REF_PREFIX="hg38-camo_mask" #Prefix given to each of the masked genomes: (will be suffixed with .ploidy_X.fa depending on ploidy)
##################################################################

bash create_camo_mask.sh \
	$ILL100_ALIGN_TO_BED \
	$REF \
	$MASK_RESULT_DIR \
	$MASK_REF_PREFIX
