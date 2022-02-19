#!/bin/bash

# Enable bash debugging to log all commands
set -x

# Receiving a string list of all input .bed files that
# were output by DRF in the previous step.
LOW_MAPQ_DIR_BED_STRING_LIST=$1

# Converting the string of input .bed files to a
# bash array.
LOW_MAPQ_DIR_BED_LIST=($LOW_MAPQ_DIR_BED_STRING_LIST)

# The prefix for output files
RESULT_PREFIX=$2

# The number of threads to split this task into
THREADS=$3

function combine() {
	SUFFIX=$1
	LOW_DEPTH_DIR_OUT=$2
	LOW_MAPQ_DIR_OUT=$3

	COMBINE_INPUT="${suffix}.input"
	find . -name "$suffix" > $COMBINE_INPUT

	touch $LOW_DEPTH_DIR_OUT
	touch $LOW_MAPQ_DIR_OUT
	combine_DRF_output.py \
		$COMBINE_INPUT \
		$LOW_DEPTH_DIR_OUT \
		$LOW_MAPQ_DIR_OUT
}

# total_lines=$(wc -l $drf_files | awk '{print $1}')

# Count the number of lines in the first .bed file to estimate
# the number of lines across all input .bed files.
total_lines=$(zcat ${LOW_MAPQ_DIR_BED_LIST[0]} | wc -l)
nline_per_batch=$(( ($total_lines + $THREADS - 1) / $THREADS ))

# For each low_mapq_bed file provided, split it into THREADS
# of $nline_per_batch and output the files into a unique directory
# for each sample.
for low_mapq_bed in ${LOW_MAPQ_DIR_BED_LIST[@]}
do
	base=$(basename $low_mapq_bed)
	sample=${base%%.*}
	sample_dir="$sample/"
	
	mkdir $sample_dir

	zcat $low_mapq_bed | split -l $((nline_per_batch)) - $sample_dir &

done
wait

LOW_DEPTH_DIR="low_depth"
LOW_MAPQ_DIR="low_mapq"

mkdir $LOW_DEPTH_DIR
mkdir $LOW_MAPQ_DIR

# All input .bed files should have all positions for the reference
# genome because we ran DRF with '--min-depth == -1' and '--min-mapq-mass == -1',
# so all samples should have the exact same number of lines in their respective
# .bed files.
for split_file in $sample_dir/*
do
	suffix=$(basename $split_file)	
	time combine $suffix "${LOW_DEPTH_DIR}/$suffix.gz" "${LOW_MAPQ_DIR}/$suffix.gz" &
done
wait

LOW_DEPTH_DIR_OUT=${RESULT_PREFIX}.combined.dark.low_depth.bed.gz
LOW_MAPQ_DIR_OUT=${RESULT_PREFIX}.combined.dark.low_mapq.bed.gz

echo -e "chrom\tstart\tend\tavg_nMapQBelowThreshold\tavg_depth\tavg_percMapQBelowThreshold" | gzip > $LOW_DEPTH_DIR_OUT
echo -e "chrom\tstart\tend\tavg_nMapQBelowThreshold\tavg_depth\tavg_percMapQBelowThreshold" | gzip > $LOW_MAPQ_DIR_OUT

# The files were written via gzip, so can simply be catted together.
cat ${LOW_DEPTH_DIR}/*  >> $LOW_DEPTH_DIR_OUT
cat ${LOW_MAPQ_DIR}/*   >> $LOW_MAPQ_DIR_OUT

# TODO: Add global setting to determine whether to cleanup scratch.
# rm -rvf $TMP_DIR
