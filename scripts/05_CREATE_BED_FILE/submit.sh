#!/bin/bash 

## 05_CREATE_BED_FILE
## Step 5 of our analysis, contains scripts for making the bedfiles and tables
## described in our manuscript. It takes the low_depth and low_mapq bed files
## created in Step 2 and intersects them with the annotation bed created in step
## 4 and quantifies what percentage of genome wide / gene bodies / CDS bases are 
## dark, stratifying by both GENCODE biotype and Gene Body element 

## This script also calculate camouflaged regions by blatting the low_mapq gene body
## regions and mapping them to each other. Low mapq gene body regions that have at
## least 98% sequence identity to another low mapq gene body region are considered camo

## the annotation file contain intersect of different dark/camo beds with gene annotations.
## the percent_dark_gene file contians list of genes and what percentage dark they are
## the percent_dark_biotypes file contains list of GENCODE biotypes and how many dark bases for each biotype
## the percent_dark_coding_region file contains list of protein_coding gene regions (UTRs, CDS,
## introns, etc.) and how many dark bases there are corresponding to each region

## Output from this step is used to run our camouflage gene short read pipeline to recover
## variants. The camo.align_to file is used to mask the human reference, the camo.realign bed file
## is used to realign reads from the original bam, and the camo.gatk file is passed into 
## GATK to define the interval where variants can possibly be called

###################################################
## TODO: Fill in Correct Paths                   ##
###################################################
REF="Homo_sapiens.GRCh38.fa" # hg38 reference
ANNO="../results/annotation/hg38/Homo_sapiens.GRCh38.93.annotation.bed" # hg38 gene annotation bed output from Step 4

ILL100_LOW_COV="../results/illuminaRL100/combined/IlluminaRL100.hg38.combined.dark.low_depth.bed" #Combined ill100 depth output from step 2
ILL100_LOW_MAPQ="../results/illuminaRL100/combined/IlluminaRL100.hg38.combined.dark.mapq.bed" #Combined ill100 mapq output from step 2
ILL100_RESULT_DIR="../results/hg38/illuminaRL100"

ILL250_LOW_COV="../results/illuminaRL250/combined/IlluminaRL250.hg38.combined.dark.low_depth.bed" #Combined ill250 depth output from step 2
ILL250_LOW_MAPQ="../results/illuminaRL250/combined/IlluminaRL250.hg38.combined.dark.mapq.bed" #Combined ill250 mapq output from step 2
ILL250_RESULT_DIR="../results/hg38/illuminaRL250"

TENX_LOW_COV="../results/10X/10X.hg38.dark.low_depth.bed" # 10X output from step 2
TENX_LOW_MAPQ="../results/10X/10X.hg38.dark.low_mapq.bed"
TENX_RESULT_DIR="../results/hg38/10X"

ONT_LOW_COV="../results/ONT/ONT.hg38.dark.low_depth.bed" #ONT output from step 2
ONT_LOW_MAPQ="../results/ONT/ONT.hg38.dark.low_mapq.bed"
ONT_RESULT_DIR="../results/hg38/ONT"
PB_LOW_COV="../results/PacBio/PacBio.hg38.dark.low_depth.bed" # PacBio output from step 2 PB_LOW_MAPQ="../results/PacBio/PacBio.hg38.dark.low_mapq.bed"
PB_RESULT_DIR="../results/hg38/PacBio"
#########################################################

##Illumina RL100
qsub camo_gene_pipeline.ogs \
	-d $ILL100_LOW_COV \
	-m $ILL100_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s illuminaRL100 \
	-v hg38 \
	-t 16 \
	-r $ILL100_RESULT_DIR

##Illumina RL250
qsub camo_gene_pipeline.ogs \
	-d $ILL250_LOW_COV \
	-m $ILL250_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s illuminaRL250 \
	-v hg38 \
	-t 16 \
	-r $ILL250_RESULT_DIR

##PacBio
qsub camo_gene_pipeline.ogs \
	-d $PB_LOW_COV \
	-m $PB_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s pacbio \
	-v hg38 \
	-t 16 \
	-r $PB_RESULT_DIR

## ONT
qsub camo_gene_pipeline.ogs \
   -d $ONT_LOW_COV \
   -m $ONT_LOW_MAPQ \
   -g $REF \
   -a $ANNO \
   -s ont \
   -v hg38 \
   -t 16 \
   -r $ONT_RESULT_DIR

## 10X
qsub camo_gene_pipeline.ogs \
	-d $TENX_LOW_COV \
	-m $TENX_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s 10x \
	-v hg38 \
	-t 16 \
	-r $TENX_RESULT_DIR
