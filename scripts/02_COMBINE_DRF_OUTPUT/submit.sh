#!/bin/bash

## 02_COMBINE_DRF_OUTPUT

## Step 2 of our analysis is to average together the 10 DRF
## outputs for each of the illumina ReadLength data. The python script
## combine_DRF_output.py does this as well as seperating the DRF output into
## a low_mapq file where MAPQ mass is ≥ 90% and a low_depth file where 
## depth is ≤ 5

## Another script split_DRF_output.ogs is used on the long- and linked-read data
## sets to seperate those DRF outputs into a low_mapq and low_depth bed

#################################
##TODO: Fill in Correct Paths  ##
#################################
#Dir that contains DRF output files from the Ill100 run in Step 1
#Dir should ONLY contain this DRF output and all files should have suffix *.low_depth.dark.bed
ILL100_DRF_OUTDIR="../results/illuminaRL100/DRF/hg38/" 
ILL100_RESULT_DIR="../results/illuminaRL100/combined/hg38/" #Where to put write combined output
ILL100_PREFIX="IlluminaRL100.hg38.combined" #Output file prefix

ILL250_DRF_OUTDIR="../results/illuminaRL250/DRF/hg38/"
ILL250_RESULT_DIR="../results/illuminaRL250/combined/hg38/"
ILL250_PREFIX="IlluminaRL250.hg38.combined"

PB_DRF_OUTDIR="../results/PacBio/DRF/hg38" #PacBio DRF output file from step 1
PB_RESULT_DIR="../results/PacBio"
PB_PREFIX="PacBio.hg38"

TENX_DRF_OUTDIR="../results/illuminaRL250/DRF" #10X DRF output file from step 1
TENX_RESULT_DIR="../results/10X"
TENX_PREFIX="10X.hg38"

ONT_DRF_OUT="../results/ONT/DRF/hg38" #ONT DRF output file from step 1
ONT_RESULT_DIR="../results/ONT"
ONT_PREFIX="ONT.hg38"
##################################

#combine Illumina RL100 Samples
bash combine_DRF.ogs $ILL100_DRF_OUTDIR $ILL100_RESULT_DIR $ILL100_PREFIX

#combine Illumina RL250 samples
bash combine_DRF.ogs $ILL250_DRF_OUTDIR $ILL250_RESULT_DIR $ILL250_PREFIX

# split PacBio output into Low Depth and Low Mapq Bed
qsub split_DRF_output.ogs $PB_DRF_OUT $PB_RESULT_DIR $PB_PREFIX

# split 10X output
qsub split_DRF_output.ogs $TENX_DRF_OUTDIR $TENX_RESULT_DIR $TENX_PREFIX

# split ONT output 
qsub split_DRF_output.ogs $ONT_DRF_OUT $ONT_RESULT_DIR $ONT_PREFIX
