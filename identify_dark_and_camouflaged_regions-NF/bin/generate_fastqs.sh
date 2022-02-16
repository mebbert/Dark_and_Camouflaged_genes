#!/bin/bash

# Enable bash debugging to log all commands
set -x


sample_input=$1
n_reads_per_run=$2
original_ref=$3
n_threads=$4
mem_per_thread=$5


echo "`date`: job running on `hostname`"
echo "`date`: Generate .fastqs for: $sample_input"


####################################
# Obtain sample name from .bam @RG #
####################################
regex="SM:([A-Za-z0-9_\-]+)"
RG=$(samtools view -H $sample_input | grep '^@RG' | tail -1)
[[ $RG =~ $regex ]] 
sample_name=${BASH_REMATCH[1]}



#####################################################
# Test whether n_reads_per_run is an even number.   #
# Must be even to make sure no read pairs are split #
# between fastq files.                              #
#####################################################
if [ $((n_reads_per_run%2)) -ne 0 ]
then
  n_reads_per_run=$((n_reads_per_run+1))
  echo "WARNING: the number of reads per run must be even. Adding one."
  echo "N reads per run now equals: $n_reads_per_run"
fi



########################################
# Get the original extension .(cr|b)am #
########################################
filename=$(basename -- "$sample_input")
extension="${filename##*.}"


# Replace actual tab characters with string representation of tab characters ('\t')
RG=${RG//	/\\\\t}




#####################
# Sort .bam by name #
#####################

# Will use .bam for intermediate files regardless of the specified output
# format (for simplicity).
sample_sorted_by_name="${sample_name}.name_sorted.bam"

echo "`date`: Sample data sorted by name will be written to $sample_sorted_by_name"

# Sorting the sample file by name to extract read pairs for realignment
echo "`date`: Sorting $sample_input by name"
if ! time samtools sort -@ $n_threads -n -m ${mem_per_thread}G $sample_input > $sample_sorted_by_name; then
	echo "ERROR: sort failed"
	exit 1
fi



########################################
# Extract reads to paired .fastq files #
########################################

# Defining the fastq names
# fq1="${sample_name}_R1.fastq"
# fq2="${sample_name}_R2.fastq"
fq="${sample_name}.interleaved_R1_R2.fastq"

# Export CRAM_REFERENCE so bamtofastq has it when the input is a .cram file
# export CRAM_REFERENCE=$original_ref

# Extracting all reads from the bam/cram and converting to fastq for
# realignment
# echo "`date`: Converting to Fastq:"
# if ! time bedtools bamtofastq -i $sample_sorted_by_name -fq $fq1 -fq2 $fq2; then
# 	echo "ERROR: bamtofastq failed"
# 	exit 1
# fi


#############################
# Create interleaved .fastq #
#############################
echo "`date`: Converting to Fastq:"
if ! time samtools fastq -o $fq -s /dev/null -0 /dev/null --reference $original_ref $sample_sorted_by_name; then
	echo "ERROR: bamtofastq failed"
	exit 1
fi

# Fix bad formatting this sample's read names (uses old format that makes bwa choke). This was used
# for one specific sample from dbGap that was originally aligned many years ago and had some funky
# read names that bwa choked on. I'm leaving this here just in case it's ever useful. This removes
# the '/1' and '/2' entirely.
#
# The variable names above and below ($fq1 and $fq2) would need to be modified to plug this back in
# to the script. Otherwise will overwrite $fq1 and $fq2 with nothing.
# awk '{ if(NR%4==1) print gensub(/(^@[[:digit:]]*)(\/[[:digit:]]$)/, "\\1", "g", $0); else print $1}' $tmp_fq1 > $fq1
# awk '{ if(NR%4==1) print gensub(/(^@[[:digit:]]*)(\/[[:digit:]]$)/, "\\1", "g", $0); else print $1}' $tmp_fq2 > $fq2


################
# Split fastqs #
################

# Split each .fastq into sets of '$n_reads_per_run'. There are 4 lines for every
# read in a .fastq file, so multiply '$n_reads_per_run' by 4.
split_dir="split_fqs"
mkdir -p $split_dir

n_lines_per_read=4
n_lines=$(( n_reads_per_run * n_lines_per_read ))

if ! time split --additional-suffix=".fastq" -dl $n_lines $fq "${split_dir}/${sample_name}.interleaved_R1_R2.split."; then
	echo "ERROR: split failed"
	exit 1
fi

# if ! time split --additional-suffix=".fastq" -dl $n_lines $fq2 "${split_dir}/${sample_name}_R2.split."; then
# 	echo "ERROR: split failed"
# 	exit 1
# fi


###########################################################
# Store sample name, .fastq, and RG to file for alignment #
###########################################################
tuple_file="${sample_name}.tuples.txt"

# Create empty file
> $tuple_file
for fq_file in $PWD/${split_dir}/*
do
	echo "${sample_name},$fq_file,$RG" >> $tuple_file
done




echo "`date` DONE!"
