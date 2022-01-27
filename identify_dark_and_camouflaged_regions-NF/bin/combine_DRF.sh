#!/bin/bash

# Enable bash debugging to log all commands
set -x

# Receiving a string list of all input .bed files that
# were output by DRF in the previous step.
LOW_MAPQ_BED_STRING_LIST=$1

# Converting the string of input .bed files to a
# bash array.
LOW_MAPQ_BED_LIST=($LOW_MAPQ_BED_STRING_LIST)

# The prefix for output files
RESULT_PREFIX=$2

# TMP_DIR="tmp/${JOB_ID}"
# mkdir -p $TMP_DIR

function combine() {
	SUFFIX=$1
	LOW_DEPTH_OUT=$2
	LOW_MAPQ_OUT=$3

	# COMBINE_INPUT="${TMP_DIR}/${suffix}.input"
	COMBINE_INPUT="${suffix}.input"
	# find $TMP_DIR -name "$suffix" > $COMBINE_INPUT
	find . -name "$suffix" > $COMBINE_INPUT

	touch $LOW_DEPTH_OUT
	touch $LOW_MAPQ_OUT
	combine_DRF_output.py \
		$COMBINE_INPUT \
		$LOW_DEPTH_OUT \
		$LOW_MAPQ_OUT
}

# total_lines=$(wc -l $drf_files | awk '{print $1}')

# Count the number of lines in the first .bed file to estimate
# the number of lines across all input .bed files.
total_lines=$(wc -l ${LOW_MAPQ_BED_LIST[0]} | awk '{print $1}')
chunks=16
nline=$(( ($total_lines + $chunks - 1) / $chunks ))

# For each low_mapq_bed file provided, split it into chunks
# of $nline and output the files into a unique directory
# for each sample.
for low_mapq_bed in ${LOW_MAPQ_BED_LIST[@]}
do
	base=$(basename $low_mapq_bed)
	sample=${base%%.*}
	# sample_dir="${TMP_DIR}/$sample/"
	sample_dir="$sample/"
	
	mkdir $sample_dir

	split -l $((nline)) $low_mapq_bed $sample_dir &
done
wait

# LOW_DEPTH="${TMP_DIR}/low_depth"
# LOW_MAPQ="${TMP_DIR}/low_mapq"

LOW_DEPTH="low_depth"
LOW_MAPQ="low_mapq"

mkdir $LOW_DEPTH
mkdir $LOW_MAPQ

# All input .bed files should have all positions for the reference
# genome because we ran DRF with '--min-depth == -1' and '--min-mapq-mass == -1',
# so all samples should have the exact same number of lines in their respective
# .bed files.
for split_file in $sample_dir/*
do
	suffix=$(basename $split_file)	
	time combine $suffix "${LOW_DEPTH}/$suffix" "${LOW_MAPQ}/$suffix" &
done
wait

LOW_DEPTH_OUT=${RESULT_PREFIX}.combined.dark.low_depth.bed
LOW_MAPQ_OUT=${RESULT_PREFIX}.combined.dark.low_mapq.bed

echo -e "chrom\tstart\tend\tavg_nMapQBelowThreshold\tavg_depth\tavg_percMapQBelowThreshold" > $LOW_DEPTH_OUT
echo -e "chrom\tstart\tend\tavg_nMapQBelowThreshold\tavg_depth\tavg_percMapQBelowThreshold" > $LOW_MAPQ_OUT
cat ${LOW_DEPTH}/* >> $LOW_DEPTH_OUT
cat ${LOW_MAPQ}/* >> $LOW_MAPQ_OUT

# TODO: Add global setting to determine whether to cleanup scratch.
# rm -rvf $TMP_DIR
