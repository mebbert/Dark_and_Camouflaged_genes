#!/bin/bash
#$ -cwd  ## start job in current working directory (where it's submitted from)
#$ -N combineVariants_ADSP  ## job name
#$ -q <queue>  ## queue
#$ -M email@institution.edu  ## email
#$ -pe threaded 1
#$ -l h_vmem=40G  ## Memory request
#$ -notify  ## tells OGS to allow a 'trap'
#$ -j y

# source .bash_profile because OGS is lame
source $HOME/.bash_profile

VCF_DIR=$1
REF=$2
GATK_JAR=$3
PICARD_JAR=$4
RESULT_DIR=$5
GATK="java -Xmx12G -Xms12G -cp $GATK_JAR"

diff_ploidy_input=""
for dir in $VCF_DIR/*
do
	ploidy=$(basename $dir)
	input_vcfs=$(find $dir -name "*.vcf" | awk '{print "-V "$1}')
	out="${RESULT_DIR}/full_cohort.ADSP.camo_genes.genotyped.hg38.${ploidy}.vcf"
	diff_ploidy_input="${diff_ploidy_input} I=${out}"
	$GATK org.broadinstitute.gatk.tools.CatVariants \
		-R $REF \
		-out $out \
		$input_vcfs
done

PICARD="java -Xmx32G -Xms32G -jar $PICARD_JAR"
$PICARD SortVcf\
	$diff_ploidy_input \
	O=${RESULT_DIR}/full_cohort.ADSP_WES.camo_genes.genotyped.hg38.combined.vcf
