#!/bin/bash

###################################
##TODO: fill in correct paths    ##
###################################
ref="Homo_sapiens.GRCh37.dna.primary_assembly.sorted.fa" #b37 reference karotypically sorted
DRF_jar="DarkRegionFinder.jar" #Path to DRF jar

ILL100_BAM_LIST="ADSP_rl100_new_bam.list" #List of paths to the 10 ADSP 100bp readlength bams, one path per line
ILL100_RESULT_DIR="../../results/illuminaRL100/DRF" #where to write DRF output for the 10 illRl100 samples

ILL250_BAM_LIST="bam_list.b37.250.txt" #List of paths to the 10 1kGenome 250bp readlength bams
ILL250_RESULT_DIR="../../results/illuminaRL250/DRF" 

PB_BAM="HG005_PacBio_GRCh37.bam" #path to PacBio bam 
PB_RESULT_DIR="../../results/PacBio/DRF"

TENX_BAM="HG00512_possorted_bam_10x.bam" #path to 10x bam
TENX_RESULT_DIR="../../results/10X/DRF"

ONT_BAM="cliveome2_plus_HG002.bam" # path to ONT bam
ONT_RESULT_DIR="../../results/ONT/DRF"
###################################

#submit jobs for illumina read length 100
while read bam
do
	qsub run_DRF.ogs $bam $ref $DRF_jar $ILL100_RESULT_DIR
done < $ILL100_BAM_LIST

#submit jobs for illumina read length 250
while read bam
do
	qsub run_DRF.ogs $bam $ref $DRF_jar $ILL250_RESULT_DIR
done < $ILL250_BAM_LIST

#Submit PacBio job
qsub run_DRF.ogs $PB_BAM $ref $DRF_jar $PB_RESULT_DIR

#Submit 10x job
qsub run_DRF.ogs $TENX_BAM $ref $DRF_jar $TENX_RESULT_DIR

#Submit ONT job
qsub run_DRF.ogs $ONT_BAM $ref $DRF_jar $ONT_RESULT_DIR
