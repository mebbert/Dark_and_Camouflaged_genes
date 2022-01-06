#!/bin/bash


## 01_RUN_DRF
## Step 1 of our analysis was running The DarkegionFinder Algorithm on the hg38 bams 
## aligned in step 0. In short, the DRF algorithm walks across the genome and for each
## counts how many reads align and what percentage of those reads are elow a MAPQ value of 10

## DRF is ran on all sequencing technologies. Parameters are set so every single position
## of the genome is output with depth and MAPQ values. Will be used later to average together
## the 10 samples from each illumina tech and to calculate median coverage of the techs


###################################
##TODO: fill in correct paths    ##
###################################
ref="Homo_sapiens.GRCh38.fa" #hg38 no alt reference karotypically sorted
DRF_jar="DarkRegionFinder.jar" #Path to DRF jar

ILL100_BAM_LIST="RL100_bams.hg38.list" #List of paths to the 10 ADSP 100bp readlength bams, aligned in step 00
ILL100_RESULT_DIR="../results/illuminaRL100/DRF/hg38" #where to write DRF output for the 10 illRl100 samples

#ILL250_BAM_LIST="RL250_bams.hg38_no_alt.list" #List of paths to the 10 1kGenome 250bp readlength bams, aligned in step 00
#ILL250_RESULT_DIR="../results/illuminaRL250/DRF/hg38" 

#PB_BAM="HG005_PacBio_GRCh38.bam" #path to PacBio bam 
#PB_RESULT_DIR="../../results/PacBio/DRF/hg38"

#TENX_BAM="HG00512_possorted_bam.hg38.bam" #path to 10x bam
#TENX_RESULT_DIR="../../results/10X/DRF/hg38"

#ONT_BAM="cliveome3_hg38.bam" # path to ONT bam
#ONT_RESULT_DIR="../../results/ONT/DRF/hg38"
###################################

#submit jobs for illumina read length 100

while read bam
do
	sbatch submitWithSingularity.slurm $bam
done < $ILL100_BAM_LIST

#while read bam
#do
#	qsub run_DRF.ogs $bam $ref $DRF_jar $ILL100_RESULT_DIR
#done < $ILL100_BAM_LIST

#submit jobs for illumina read length 250
#while read bam
#do
#	qsub run_DRF.ogs $bam $ref $DRF_jar $ILL250_RESULT_DIR
#done < $ILL250_BAM_LIST

##Submit PacBio job
#qsub run_DRF.ogs $PB_BAM $ref $DRF_jar $PB_RESULT_DIR

##Submit 10x job
#qsub run_DRF.ogs $TENX_BAM $ref $DRF_jar $TENX_RESULT_DIR

##Submit ONT job
#qsub run_DRF.ogs $ONT_BAM $ref $DRF_jar $ONT_RESULT_DIR
