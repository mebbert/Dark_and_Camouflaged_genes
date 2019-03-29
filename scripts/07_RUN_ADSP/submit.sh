#!/bin/bash


## Step 07_RUN_ADSP

## In step 07 of our analysis, we show script on how to rescue variants
## from camouflaged regions. First, we take the realign file from step 5
## which describes where all camo regions are located. And we extract low MAPQ reads
## from these regions from original bams. We realign these reads to the masked genomes
## created in step 6. After realigning we call variants in the GATK interval created
## in step 5.

## We do this seperately for each repeat number, adjusting the ploidy of GATK Haplotype Caller
## for the specific region we're calling (i.e. regions with repeat number 2 we use ploidy of 4)
## After which we combineGVCFs from all samples from the same ploidy to get combined GVCFs for 
## each ploidy tested.

## We parallelize this workflow across multiple compute nodes, having each node analyze 50 bams, in 
## later steps we will combine all intermediate gVCFs into a single gVCF for the whole cohort


##################################################################
## TODO: Fill in Correct Paths                                  ##
##################################################################
REALIGN="hg19.camo.realign.bed" # realign camo bed created in step 5 (added "chr" to all b37 camo regions so positions will be in hg19)
GATK_BED="../results/hg38/illuminaRL100/illuminaRL100.hg38.camo.GATK.bed" # ill100 GATK bed created in step 5

GATK_JAR="GenomeAnalysisTK.jar" #path to GATK jar (version 3.8)

# REF_PREFIX is equal to "${MASK_RESULT_DIR}/${MASK_REF_PREFIX}" from Step 6
REF_PREFIX="../results/hg38_camo_mask/hg38-camo_mask" # Path prefix to the camo-masked seperate ploidy references
FILTERED_BAM_LIST="ADSP.WES.bams.list" # List of paths to ADSP (or other dataset) bams that passed QC, one per line
DATA_DEST_DIR="../results/ADSP_WES" #Where to store the intermediate gVCFs results 

##################################################################

##Run Filtered Samples
numSamples=$(wc -l $FILTERED_BAM_LIST | awk '{print $1}')

# Parallelize job across compute node, each node  analyzing 50 bams
for i in $(seq 1 50 $numSamples)
do
	qsub run_camo_genes.ogs \
		$i,$((i + 49)) \
		$REALIGN \
		$GATK_BED \
		$GATK_JAR \
		$REF_PREFIX \
		$FILTERED_BAM_LIST \
		$DATA_DEST_DIR 
done
