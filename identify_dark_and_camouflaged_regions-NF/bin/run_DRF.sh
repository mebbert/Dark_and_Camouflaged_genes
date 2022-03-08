#!/bin/bash

# Enable bash debugging to log all commands
set -x

sample_name=$1
sample_input_file=$2 #Bam input into DRF program
align_to_ref=$3 # Human Reference genome used in the alignment 
DRF_jar=$4 #Path to DRF jar
intervals=$5

#Regex to Extract Sample name from Bam file
# sm_regex="SM:([A-Za-z0-9_\-]+)"
# RG=$(samtools view -H $sample_input_file | grep '^@RG' | head -1)
# [[ $RG =~ $sm_regex ]]
# name=${BASH_REMATCH[1]}

sample_dir="${sample_name}_mini_beds"

mkdir -p $sample_dir

echo "`date`: Running sample $sample_name"

## Calculate depth and mapq mass for every position in the genome
## Allows us to calculate average depth and average across samples
min_depth=-1 #MAKE MIN_DEPTH = -1 ensures no bases are printed in LOW_COV_BED
min_mapq_mass=-1 #MAKE MAPQ MASS = -1 ensures EVERY POSITION GETS PRINTED IN LOW_MAPQ_BED
mapq_thresh=9 #gatk recommended cutoff
low_cov_bed="/dev/null" # Nothing output to low_cov_bed since min_depth = -1
low_mapq_bed="${sample_dir}/${sample_name}.min_depth_${min_depth}.min_mapq_mass_${min_mapq_mass}.mapq_thresh_${mapq_thresh}.dark.low_mapq.bed"
inc_bed="/dev/null" # No need to store incomplete bases for each sample

if ! java -Xmx32G -jar $DRF_jar \
		-i "${sample_input_file}" \
		--human-ref "${align_to_ref}" \
		--min-region-size 1 \
		--mapq-threshold "${mapq_thresh}" \
		--min-mapq-mass "${min_mapq_mass}" \
		--min-depth "${min_depth}" \
		--region-exclusivity \
		--low-coverage-bed-output "${low_cov_bed}" \
		--low-mapq-bed-output "${low_mapq_bed}" \
		--incomplete-bed-output "${inc_bed}" \
		--interval-list "${intervals}"; then

	echo "ERROR: DRF failed. Check logs for details."
	exit 1
fi
