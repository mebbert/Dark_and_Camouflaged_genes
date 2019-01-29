#!/bin/bash

############################################################
## TODO: Fill in correct paths                            ##
############################################################
VCF_DATA_DIR="../../results/ADSP/genotyped_by_region" #Result directory from step 8 where intermediate VCFs are written
REF="Homo_sapiens.GRCh37.dna.primary_assembly.sorted.fa" #b37 reference
GATK_JAR="GenomeAnalysisTK.jar" #Path to GATK jar
CAMO_BED="../../results/IlluminaRL100/illuminaRL100.b37.camo.bed" #ILL100 merged camo bed created in step 5
ANNO_BED="../../results/annotations/b37/Homo_sapiens.GRCh37.87.annotation.bed" #Gene annotation b37 created in step 4
RESULT_DIR="../../results/ADSP" #Where to write results
############################################################

bash catVariants.sh \
	$VCF_DATA_DIR \
	$REF \
	$GATK_JAR \
	$RESULT_DIR

INPUT_VCF="${RESULT_DIR}/full_cohort.ADSP.camo_genes.genotyped.b37.combined.vcf"
ANNOTATED="${RESULT_DIR}/annotated.vcf"
python calc_genotype_annotations.py $INPUT_VCF > $OUT

INTERSECTED="${RESULT_DIR}/annotated.intersected.CDS.bed"
bash intersect.sh \
	$ANNOTATED \
	$CAMO_BED \
	$ANNO_BED \
	> $INTERSECTED

FILTERED="${RESULT_DIR}/annotated.intersected.CDS.filtered.bed"
python filter.py $INTERSECTED > $FILTERED

VARIANT_METRICS="${RESULT_DIR}/annotated.intersected.CDS.filtered.variant_metrics.txt"
python extractVariantQualityMetrics.py $FILTERED > $VARIANT_METRICS
