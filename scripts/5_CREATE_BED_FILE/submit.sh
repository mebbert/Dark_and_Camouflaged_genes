#!/bin/bash 

###################################################
## TODO: Fill in Correct Paths                   ##
###################################################
REF="Homo_sapiens.GRCh37.dna.primary_assembly.sorted.fa" # b37 reference, Karyotypically sorted
ANNO="../../results/annotations/b37/Homo_sapiens.GRCh37.87.annotation.bed" # Gene annotation bed output from Step 4

ILL100_LOW_COV="../../results/illuminaRL100/IlluminaRL100.b37.combined.dark.low_depth.bed" #Combined ill100 DRF output from step 2
ILL100_LOW_MAPQ="../../results/illuminaRL100/IlluminaRL100.b37.combined.dark.low_mapq.bed"
Ill100_RESULT_DIR="../../results/illuminaRL100/b37"

ILL250_LOW_COV="../../results/illuminaRL250/illuminaRL250.b37.combined.dark.low_depth.bed" #Combined Illumina RL250 DRF output from step 2
ILL250_LOW_MAPQ="../../results/illuminaRL250/illuminaRL250.b37.combined.dark.low_mapq.bed"
ILL250_RESULT_DIR="../../results/illuminaRL250/b37"

TENX_LOW_COV="../../results/10X/10X.b37.dark.low_depth.bed" # 10X output from step 2
TENX_LOW_MAPQ="../../results/10X/10X.b37.low_mapq.bed"
TENX_RESULT_DIR="../../results/10X/b37"

ONT_LOW_COV="../../results/ONT/ONT.b37.dark.low_depth.bed" #ONT output from step 2
ONT_LOW_MAPQ="../../results/ONT/ONT.b37.dark.low_mapq.bed"
ONT_RESULT_DIR="../../results/ONT/b37"

PB_LOW_COV="../../results/PacBio/PacBio.b37.dark.low_depth.bed" # PacBio output from step 2
PB_LOW_MAPQ="../../results/PacBio/PacBio.b37.dark.low_mapq.bed"
PB_RESULT_DIR="../../results/PacBio/b37"
#########################################################

#Illumina RL100
qsub camo_gene_pipeline.sh \
	-d $ILL100_LOW_COV \
	-m $ILL100_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s illuminaRL100 \
	-v b37 \
	-t 16 \
	-r $ILL100_RESULT_DIR

#Illumina RL250
qsub camo_gene_pipeline.sh \
	-d $ILL250_LOW_COV \
	-m $ILL250_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s illuminaRL250 \
	-v b37 \
	-t 16 \
	-r $ILL1250_RESULT_DIR

#PacBio
qsub camo_gene_pipeline.sh \
	-d $PB_LOW_COV \
	-m $PB_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s pacbio \
	-v b37 \
	-t 16 \
	-r $PB_RESULT_DIR

# ONT
qsub camo_gene_pipeline.sh \
	-d $ONT_LOW_COV \
	-m $ONT_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s ont \
	-v b37 \
	-t 16 \
	-r $ONT_RESULT_DIR

# 10X
qsub camo_gene_pipeline.sh \
	-d $TENX_LOW_COV \
	-m $TENX_LOW_MAPQ \
	-g $REF \
	-a $ANNO \
	-s 10x \
	-v b37 \
	-t 16 \
	-r $TENX_RESULT_DIR
