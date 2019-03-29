#!/bin/bash

## 09_FIND_FALSE_POSITIVES

## in step 9 of our analysis, we create a list of all possible reference based 
## artifacts. To do this we take all the camo CDS regions with repeat number ≤ 5, 
## and blat them against the whole genome. 

## Any DNA sequence from a hit with sequence identity ≥ 98% were locally
## aligned back to the query seqeunce using Bio.pairwise2.
## Any mismatches or gaps in the aligned sequence were converted into variant
## positions and output to a bed file that lists positions of reference-based-artifacts

## reference-based artifcats are false positive variant calls created when 
## two regions in a camo set are not completely identical, so when a read from one
## region is forced to align to the other, it will create an artificial variant

## In a later filtering any variant called by GATK from one of these positions is 
## filtered out

###########################################################
### Fill In Correct Paths:                              ###
###########################################################
CAMO_ANNOTATION="../results/hg38_no_alt/illuminaRL100/illuminaRL100.hg38_no_alt.camo_annotations.txt"
ALIGN_TO="../results/hg38_no_alt/illuminaRL100/illuminaRL100.hg38_no_alt.camo.align_to.sorted.bed"
REF="/research/labs/moleneurosci/petrucelli/m199515/NGS_resources/references/refdata-GRCh38-2.1.0/fasta/genome.fa"

RESULT="../results/ADSP_WES/reference_based_artifacts.hg38_no_alt.bed" #where to write reference-based artifacts bed
###########################################################

qsub get_false_positives.ogs \
	$CAMO_ANNOTATION \
	$ALIGN_TO \
	$REF \
	$RESULT
