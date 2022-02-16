#!/bin/bash

# Enable bash debugging to log all commands
set -x

sample_alignment_file=$1
sample_name=$2
out_format_param=$3
n_threads=$4
mem_per_thread=$5


echo "`date`: job running on `hostname`"


# Determine the requested output file format (.bam or .cram). The ',,' notation
# converts the string to lowercase.
out_format=""
if [[ ${out_format_param,,} == *"bam"* ]]; then
	out_format="bam"
elif [[ ${out_format_param,,} == *"cram"* ]]; then
	out_format="cram"
else
	echo "ERROR: Invalid output format specified: $out_format_param. Must be '.bam' or '.cram' (with or without the '.' character)."
	exit 1
fi



regex="([A-Za-z0-9_\-]+).unsorted.mini.(split.[0-9]+).(cram|bam|CRAM|BAM|Cram|Bam)"
[[ ${sample_alignment_file} =~ $regex ]] 
sample_name_from_file_name=${BASH_REMATCH[1]}
split_num=${BASH_REMATCH[2]}


#############################################
# Test that sample name provided and sample #
# name from file name are identical.        #
#############################################

if [[ ${sample_name} != ${sample_name_from_file_name} ]]; then
	echo "ERROR: Sample name provided to script does not match" \
		  " sample name obtained from input file." \
		  "Provided name: $sample_name" \
		  "Name from file: $sample_name_from_file_name" \
		  "File name obtained from: $sample_alignment_file"
	exit 1
fi


echo "`date`: Sorting and indexing $sample_alignment_file" 

final_mini_output_file="${sample_name}.sorted.mini.${split_num}.${out_format}"



###################################
# Sorting the aligned sample file #
###################################

echo "`date` Sorting final sample $final_mini_output_file"
if [ "$out_format" = "bam" ]; then
	if ! time samtools sort -@ $n_threads -m ${mem_per_thread}G $sample_alignment_file > $final_mini_output_file; then
		echo "ERROR: Sorting final sample $out_format failed"
		exit 1
	fi
else
	if ! time samtools sort -O cram -@ $n_threads -m ${mem_per_thread}G $sample_alignment_file > $final_mini_output_file; then
		echo "ERROR: Sorting final sample $out_format failed"
		exit 1
	fi
fi



#################################
# Indexing the sorted .bam/cram #
#################################
echo "`date` Indexing final $out_format file ($final_mini_output_file)"
if ! time samtools index $final_mini_output_file; then
	echo "ERROR: Indexing final $out_format failed"
	exit 1
fi


# Validating the final .bam/cram for output.
echo "`date` Validating $out_format"
if ! time samtools quickcheck $final_mini_output_file; then
	echo "quickcheck failed"
	exit 1
else
	echo "quickcheck happy!"
fi




#####################################################################
# Create tuple file containing sample name, aligned file, and index #
#####################################################################
tuple_file="${sample_name}.${split_num}.tuples.txt"

index_file=(${PWD}/${final_mini_output_file}.*)
echo "${sample_name},${PWD}/${final_mini_output_file},${index_file}" > $tuple_file

echo "`date` DONE!"
#rm -rfv $TMP_DIR
