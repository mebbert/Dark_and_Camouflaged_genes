#!/bin/bash

## 08_COMBINE_AND_GENOTYPE_ADSP

## In step 8 of our analysis, we combine the intermediate gVCFs created from 
## the scatter jobs of step 7 into a single gVCF and then genotype it to create
## a cohort level vcf for each ploidy

## We found combining all 10,993 gVCFs across all camo regions to require too much
## memory so we scatter this job as well over the camo region. So we look at each camo 
## region individually, and combine all samples and genotype them for just this single region

## In a later step we will cat the Variants from all these seperate region VCF into a combined
## cohort VCF for all ploidies and for all samples across all regions.

##################################################################
## TODO: Fill in Correct paths                                  ##
##################################################################
GVCF_RESULT_DIR="../results/ADSP_WES/camo_gvcfs" # Dir created in step 7 that contains all intermediate gvcfs
GATK_BED="../results/hg38/illuminaRL100/illuminaRL100.hg38.camo.GATK.bed" # Camo GATK created in step 5
RESULT_DIR="../results/ADSP_WES" # where to write the combined vcfs
CAMO_MASK_REF_PREFIX="../results/hg38_camo_mask/hg38-camo_mask" #prefix to all the masked genomes created in step 6
GATK_JAR="GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar" # Path to GATK jar (version 3.8)
##################################################################

bash submit_combine_and_genotype.ogs \
	$GVCF_RESULT_DIR \
	$GATK_BED \
	$RESULT_DIR \
	$CAMO_MASK_REF_PREFIX \
	$GATK_JAR
