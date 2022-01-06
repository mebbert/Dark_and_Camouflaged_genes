#!/bin/bash

## 03_Calc_Bam_Metrics

## Step 3 of our analysis, we calculate metrics for the 
## alignements of step 0, quantifying coverage and read length
## for the long read technologies

## The script calc_median_depth.ogs uses the DRF output that
## print depth for every single base to calculate median and average
## depth for all bases on the main contigs 

## The script calcLengthMetrics.ogs uses the Long-read bams to to
## calculate median and N50 readlengths for each of the long read teachs

##########################################
## TODO: Fill in correct paths          ##
##########################################
ILL100_DRF_OUTDIR="../results/illuminaRL100/DRF/hg38" #Dir with illuminaRL100 output from Step 1

#ILL250_DRF_OUTDIR="../results/illuminaRL250/DRF/hg38"

#PB_DRF_OUTDIR="../results/PacBio/DRF/hg38"
#PB_BAM="HG005_PacBio_GRCh38.bam" #path to PacBio Bam

#TENX_DRF_OUTDIR="../results/10X/DRF/hg38/"

#ONT_DRF_OUTDIR="../results/ONT/DRF/hg38/"
#ONT_BAM="cliveome_v3.hg38.downsampled.bam" #path to ONT bam

RESULT_DIR="../results/bam_metrics"
###########################################

# Calculate Coverage of Bams
mkdir -p $RESULT_DIR
OUT="${RESULT_DIR}/LR_depth_metrics.txt"
echo -e "DRF_FILE\tSEQUENCER\tMEAN_DEPTH\tMEDIAN_DEPTH" > $OUT


##Submit Long read Jobs
#qsub calc_median_depth.ogs $PB_DRF_OUTDIR "PacBio_CCS" $OUT
#qsub calc_median_depth.ogs $TENX_DRF_OUTDIR "10X" $OUT
#qsub calc_median_depth.ogs $ONT_DRF_OUTDIR "ONT" $OUT

#Calculate Length Metrics of Bams
#qsub calcLengthMetrics.ogs $PB_BAM $ONT_BAM $RESULT_DIR

## Submit Ill100 Jobs
for drf_bed in $ILL100_DRF_OUTDIR/*
do
	sbatch submitWithSingularity.slurm $drf_bed "IlluminaRL100" $OUT
done

##Submit ILL250 Jobs
#for drf_bed in $ILL250_DRF_OUTDIR/*
#do
	#qsub calc_median_depth.ogs $drf_bed "IlluminaRL250" $OUT
#done
