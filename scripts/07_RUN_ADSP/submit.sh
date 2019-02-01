#!/bin/bash


##################################################################
## TODO: Fill in Correct Paths                                  ##
##################################################################
REALIGN="../../results/illuminaRL100/b37/illuminaRL100.b37.camo.realign.sorted.bed" # realign camo bed created in step 5
ALIGN_TO="../../results/illuminaRL100/b37/illuminaRL100.b37.camo.align_to.sorted.expanded_50.bed" #align_to camo bed created in step 5 and expanded by 50 in step 6

GATK_JAR="GenomeAnalysisTK.jar" #path to GATK jar (version 3.8)

# REF_PREFIX is equal to "${MASK_RESULT_DIR}/${MASK_REF_PREFIX}" from Step 6
REF_PREFIX="../../results/b37_camo_mask/b37-camo_mask" # Path prefix to the camo-masked seperate ploidy references
FILTERED_BAM_LIST="ADSP_filtered_bams.list" # List of paths to ADSP (or other dataset) bams that passed QC, one per line
DATA_DEST_DIR="../../results/ADSP" #Where to store the intermediate gVCFs results 

##################################################################

##Run Filtered Samples
numSamples=$(wc -l $FILTERED_BAM_LIST | awk '{print $1}')

# Scatter bams with each job analyzing 25 bams
for i in $(seq 1 25 $numSamples)
do
	qsub run_camo_genes.ogs \
		$i,$((i + 1)) \
		$REALIGN \
		$ALIGN_TO \
		$GATK_JAR \
		$REF_REFIX \
		$FILTERED_BAM_LIST \
		$DATA_DEST_DIR 
done
