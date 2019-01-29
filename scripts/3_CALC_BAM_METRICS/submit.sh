#!/bin/bash

##########################################
## TODO: Fill in correct paths          ##
##########################################
ILL100_DRF_OUTDIR="../../results/illuminaRL100/DRF" #Dir with illuminaRL100 output from Step 1

ILL250_DRF_OUTDIR="../../results/illuminaRL250/DRF"

PB_DRF_OUTDIR="../../results/PacBio/DRF"
PB_BAM="HG005_PacBio_GRCh37.bam" #path to PacBio Bam

TENX_DRF_OUTDIR="../../results/10X/DRF"

ONT_DRF_OUTDIR="../../results/ONT/DRF"
ONT_BAM="cliveome2_plus_HG002.bam" #path to ONT bam

RESULT_DIR="../../results"
###########################################

# Calculate Coverage of Bams
OUT="${RESULT_DIR}/depth_metrics.txt"
echo -e "DRF_FILE\tSEQUENCER\tMEAN_DEPTH\tMEDIAN_DEPTH"

## Submit Ill100 Jobs
for drf_bed in $ILL100_DRF_OUTDIR/*
do
	qsub calc_median_depth.ogs $drf_bed "IlluminaRL100" $OUT
done

##Submit ILL250 Jobs
for drf_bed in $ILL250_DRF_OUTDIR/*
do
	qsub calc_median_depth.ogs $drf_bed "IlluminaRL250" $OUT
done

##Submit Long read Jobs
qsub calc_median_depth.ogs $PB_DRF_OUTDIR/* "PacBio" $OUT
qsub calc_median_depth.ogs $TENX_DRF_OUTDIR/* "10X" $OUT
qsub calc_median_depth.ogs $ONT_DRF_OUTDIR/* "ONT" $OUT

# Calculate Length Metrics of Bams
qsub calcLengthMetrics.ogs $PB_BAM $ONT_BAM $RESULT_DIR
