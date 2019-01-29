#!/bin/bash


##################################################################
## TODO: Fill in Correct paths                                  ##
##################################################################
GVCF_RESULT_DIR="../../results/ADSP/gvcfs" # Dir created in step 7 that contains all intermediate gvcfs
ALIGN_TO_BED="../../results/illuminaRL100/b37/illuminaRL100.b37.camo.align_to.sorted.bed" # Camo align_to bed created in step 5
RESULT_DIR="../../results/ADSP" # where to write the combined vcfs
CAMO_MASK_REF_PREFIX="../../b37_camo_mask/b37-camo_mask" #prefix to all the masked genomes created in step 6
GATK_JAR="GenomeAnalysisToolKit.3.8.1.jar" # Path to GATK jar
##################################################################

bash submit_combine_and_genotype.ogs \
	$GVCF_DIR \
	$ALIGN_TO_BED \
	$RESULT_DIR \
	$CAMO_MASK_REF_PREFIX \
	$GATK_JAR
