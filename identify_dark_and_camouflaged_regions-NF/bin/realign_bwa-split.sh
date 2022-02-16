#!/bin/bash

# Enable bash debugging to log all commands
set -x


fq=$1
# fq2=$2
sample_name=$2
sample_RG=$3
align_to_ref=$4
out_format_param=$5
n_threads=$6


echo "Sample RG received: \"$sample_RG\""

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


echo "`date`: job running on `hostname`"



#######################################
# Get sample name from $fq1 file name #
#######################################

# Nextflow .splitFastq appears to add '.[0-9]+' before the file extension
regex="([A-Za-z0-9_\-]+).interleaved_R1_R2.(split.[0-9]+).fastq"
[[ $fq =~ $regex ]] 
sample_name_from_file_name=${BASH_REMATCH[1]}
split_num=${BASH_REMATCH[2]}


echo "RG: \"$sample_RG\""

#############################################
# Test that sample name provided and sample #
# name from file name are identical.        #
#############################################

if [[ ${sample_name} != ${sample_name_from_file_name} ]]; then
	echo "ERROR: Sample name provided to script does not match" \
		  " sample name obtained from .fastq file." \
		  "Provided name: $sample_name" \
		  "Name from file: $sample_name_from_file_name" \
		  "File name obtained from: $fq"
	exit 1
fi


echo "`date`: Realigning sample $sample_name" 

# Fix bad formatting this sample's read names (uses old format that makes bwa choke). This was used
# for one specific sample from dbGap that was originally aligned many years ago and had some funky
# read names that bwa choked on. I'm leaving this here just in case it's ever useful. This removes
# the '/1' and '/2' entirely.
#
# The variable names above and below ($fq1 and $fq2) would need to be modified to plug this back in
# to the script. Otherwise will overwrite $fq1 and $fq2 with nothing.
# awk '{ if(NR%4==1) print gensub(/(^@[[:digit:]]*)(\/[[:digit:]]$)/, "\\1", "g", $0); else print $1}' $tmp_fq1 > $fq1
# awk '{ if(NR%4==1) print gensub(/(^@[[:digit:]]*)(\/[[:digit:]]$)/, "\\1", "g", $0); else print $1}' $tmp_fq2 > $fq2



################################################################
# Aligning fastqs to $align_to_ref and converting to .bam file #
################################################################

tmp_sample="${sample_name}.unsorted.mini.${split_num}.${out_format}"

if [ "$out_format" = "bam" ]; then
	echo "`date` Aligning to $sample_name to ${align_to_ref}. Output will be .bam file"
	# time bwa mem -M -R "$sample_RG" -t $n_threads $align_to_ref $fq1 $fq2 | samtools view -hb - > $tmp_sample
	time bwa mem -p -M -R "$sample_RG" -t $n_threads $align_to_ref $fq | samtools view -hb - > $tmp_sample
else
	echo "`date` Aligning to $sample_name to ${align_to_ref}. Output will be .cram file"
	#time bwa mem -M -R "$sample_RG" -t $n_threads $align_to_ref $fq1 $fq2 | samtools view -C -T $align_to_ref - > $tmp_sample
	time bwa mem -p -M -R "$sample_RG" -t $n_threads $align_to_ref $fq | samtools view -C -T $align_to_ref - > $tmp_sample
fi


# Check if the exit status for every cmd in the pipe (bwa and samtools) were successful.
# Cannot simply use 'if ! <cmd>' if piping multiple cmds together. The environment
# variable $PIPESTATUS stores the status for all cmds in the pipe.
TMPSTATUS=("${PIPESTATUS[@]}")

# Check if $TMPSTATUS contains exactly two '0' (would need to add more '0' if there were more than
# two cmds.
if [[ ${TMPSTATUS[@]} == "0 0" ]];then
    echo "SUCCESS: bwa mem and samtools view completed successfully."
else
    echo "ERROR: bwa mem or samtools view failed. Check the logs for details."
    exit 1
fi




#############################################################
# Create tuple file containing sample name and aligned file #
#############################################################
tuple_file="${sample_name}.${split_num}.tuples.txt"

# Create empty file
> $tuple_file
echo "${sample_name},${PWD}/${tmp_sample}" > $tuple_file

echo "`date` DONE!"
#rm -rfv $TMP_DIR
