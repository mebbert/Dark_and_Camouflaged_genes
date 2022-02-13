#!/bin/bash

# Enable bash debugging to log all commands
set -x

path_to_alignment_files=$1
sample_name=$2
ref_tag=$3
out_format_param=$4
n_threads=$5
mem_per_thread=$6


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


file_list=$(path_to_alignment_files*.{cram,bam})

first_file="${file_list[0]}"

regex="([A-Za-z0-9_\-]+).sorted.mini.(cram|bam|CRAM|BAM|Cram|Bam)"
[[ ${first_file} =~ $regex ]] 
sample_name_from_file_name=${BASH_REMATCH[1]}


#############################################
# Test that sample name provided and sample #
# name from file name are identical.        #
#############################################

if [[ ${sample_name} != ${sample_name_from_file_name} ]]; then
	echo "ERROR: Sample name provided to script does not match" \
		  " sample name obtained from input file." \
		  "Provided name: $sample_name" \
		  "Name from file: $sample_name_from_file_name" \
		  "File name obtained from: $first_file"
	exit 1
fi



echo "`date`: Merging and indexing files located at $path_to_alignment_files" 

final_output_file="${sample_name}.${ref_tag}.sorted.merged.final.${out_format}"


###################################
# Sorting the aligned sample file #
###################################

echo "`date` Sorting final sample $final_output_file"
if [ "$out_format" = "bam" ]; then
	if ! time samtools merge -pc -@ $n_threads -m ${mem_per_thread}G -o $final_output_file ${file_list[@]}; then
		echo "ERROR: Merging mini $out_format files failed"
		exit 1
	fi
else
	if ! time samtools merge -pc -O cram -@ $n_threads -m ${mem_per_thread}G -o $final_output_file ${file_list[@]}; then
		echo "ERROR: Merging mini $out_format files failed"
		exit 1
	fi
fi



#################################
# Indexing the sorted .bam/cram #
#################################
echo "`date` Indexing final $out_format file ($final_output_file)"
if ! time samtools index $final_output_file; then
	echo "ERROR: Indexing final $out_format failed"
	exit 1
fi


##############################################
# Validating the final .bam/cram for output. #
##############################################
echo "`date` Validating $out_format"
if ! time samtools quickcheck $final_output_file; then
	echo "quickcheck failed"
	exit 1
else
	echo "quickcheck happy!"
fi


echo "`date` DONE!"
#rm -rfv $TMP_DIR
