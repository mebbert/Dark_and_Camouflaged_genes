#!/bin/bash

## 10_VARIANT_FILTERING

## Step 10 is the final analysis step of our manuscript.
## Here we combine the genotyped_by_region VCFs from all ploidies
## created in step 8, so at the end of this step we have a single VCF
## for the whole cohort

## Once we have this vcf we annotate the vcf to add an inbreeding value
## calculated treating all possible polyploidy heterozygous genotypes combined
## (i.e. we treat the genotypes 0/0/0/1, 0/0/1/1, and 0/1/1/1 as the same). During
## annotation we also filter out variants that homozygous singletons (variants only present in one
## sample, and that sample is all the way homozygous ALT 1/1/1/1)

## After annotating the vcf, we filter the VCF based on the reference-based artifacts positions
## calculated in step 9 and on QD value. We chose the QD value of 2.0 to filter variants in order
## to optimise the TiTv of the call set and maximize number of variants called

## Finally, for this QD cutoff of 2.0, and for a more stringent cutoff of 3.0, we calcualte the number of 
## variants and the number of gene regions these variants were found in.


############################################################
## TODO: Fill in correct paths                            ##
############################################################
VCF_DATA_DIR="../results/ADSP_WES/genotyped_by_region" #Result directory from step 8 where intermediate VCFs are written
REF_BASED_ARTIFACT_BED="../results/ADSP_WES/reference_based_artifacts.hg38.bed" #False Positive result bed from step 9
REF="Homo_sapiens.GRCh38.fa" #hg38 reference
GATK_JAR="GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar" #Path to GATK jar
PICARD_JAR="picard.jar" # Path to Picard Jar
ANNO_BED="../results/annotations/hg38/Homo_sapiens.GRCh38.93.annotation.hg38.bed"
RESULT_DIR="../results/ADSP_WES" #Where to write results
############################################################

bash catVariants.sh \
	$VCF_DATA_DIR \
	$REF \
	$GATK_JAR \
	$PICARD_JAR \
	$RESULT_DIR

INPUT_VCF="${RESULT_DIR}/full_cohort.ADSP_WES.camo_genes.genotyped.hg38.combined.vcf"
ANNOTATED="${RESULT_DIR}/annotated.vcf"
#python calc_genotype_annotations.py $INPUT_VCF > $ANNOTATED

FILTERED="${RESULT_DIR}/filtered.vcf"
#python filter.py $ANNOTATED $REF_BASED_ARTIFACT_BED > $FILTERED

VARIANT_METRICS="${RESULT_DIR}/filtered.variant_metrics.txt"
#python extractVariantQualityMetrics.py $FILTERED > $VARIANT_METRICS

bash getGeneNum.sh \
	$VARIANT_METRICS \
	$ANNO_BED \
	$RESULT_DIR
