#!/bin/bash

#################################
##TODO: Fill in Correct Paths  ##
#################################
#Dir that contains DRF output files from the Ill100 run in Step 1
#Dir should ONLY contain this DRF output and all files should have suffix *.low_depth.dark.bed
ILL100_DRF_OUTDIR="../../results/illuminaRL100/DRF" 
ILL100_RESULT_DIR="../../results/illuminaRL100" #Where to put write combined output
ILL100_PREFIX="IlluminaRL100.b37.combined" #Output file prefix

ILL250_DRF_OUTDIR="../../results/illuminaRL250/DRF"
ILL250_RESULT_DIR="../../results/illuminaRL250"
ILL250_PREFIX="IlluminaRL250.b37.combined"

PB_DRF_OUTDIR="../../results/PacBio/DRF" #PacBio DRF output file from step 1
PB_RESULT_DIR="../../results/PacBio"
PB_PREFIX="PacBio.b37"

TENX_DRF_OUTDIR="../../results/illuminaRL250/DRF" #10X DRF output file from step 1
TENX_RESULT_DIR="../../results/10X"
TENX_PREFIX="10X.b37"

ONT_DRF_OUTDIR="../../results/ONT/DRF" #ONT DRF output file from step 1
ONT_RESULT_DIR="../../results/ONT"
ONT_PREFIX="ONT.b37"
##################################

#combine Illumina RL100 Samples
qsub combine_DRF.ogs $ILL100_DRF_OUTDIR $ILL100_RESULT_DIR $ILL100_PREFIX

#combine Illumina RL250 samples
qsub combine_DRF.ogs $ILL250_DRF_OUTDIR $ILL250_RESULT_DIR $ILL250_PREFIX

# split PacBio output into Low Depth and Low Mapq Bed
qsub split_DRF_output.ogs $PB_DRF_OUTDIR $PB_RESULT_DIR $PB_PREFIX

# split 10X output
qsub split_DRF_output.ogs $TENX_DRF_OUTDIR $TENX_RESULT_DIR $TENX_PREFIX

# split ONT output 
qsub split_DRF_output.ogs $ONT_DRF_OUTDIR $ONT_RESULT_DIR $ONT_PREFIX
