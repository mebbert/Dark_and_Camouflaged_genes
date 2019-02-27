#!/bin/bash

## Here we take a reference genome and the Camo align_to file created
# In step 6 and mask the genome so only the align to bed regions are unmasked
# We only mask camo regions with repeat number ≤ 5, and have a seperate masked
# Genome for each ploidy we test

##################################################################
## TODO: Fill in correct paths                                  ##
##################################################################
ILL100_ALIGN_TO_BED="../../results/illuminaRL100/b37/illuminaRL100.b37.camo.align_to.sorted.bed" #Align_to camo bed created in step 5

REF="Homo_saiens.GRCh37.dna.primary_assembly.sorted.fa" #b37 reference
MASK_RESULT_DIR="../../results/b37_camo_mask" #where to write masked genomes
MASK_REF_PREFIX="b37-camo_mask" #Prefix given to each of the masked genomes: (will be suffixed with .ploidy_X.fa depending on ploidy)
##################################################################

bash create_camo_mask.sh \
	$ILL100_ALIGN_TO_BED \
	$REF \
	$MASK_RESULT_DIR \
	$MASK_REF_PREFIX
