#!/bin/bash

# Enable bash debugging to log all commands
set -x


sample_input=$1
align_to_ref=$2
ref_tag=$3
original_ref=$4
out_format_param=$5
n_threads=$6
mem_per_thread=$7


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
echo "`date`: Realign $1"

regex="SM:([A-Za-z0-9_\-]+)"
RG=$(samtools view -H $sample_input | grep '^@RG' | tail -1)
[[ $RG =~ $regex ]] 
sampleName=${BASH_REMATCH[1]}

# Replace actual tab characters with string representation of tab characters ('\t')
RG=${RG//	/\\t}

# There was one specific sample from dbGap that had a hidden 'End of Medium' character (hex 19) in
# the RG that was making the bwa command (below) choke. Drove me nuts trying to figure out why it
# kept failing. The character seemed to originate from a single quote character in "Alzheimer's".
# This command just wipes out the entire 'DS:' and replaces it with 'ADSP'. I'm leaving this here
# just in case it's ever useful.
# RG=$(echo $RG | sed -r "s/DS:Alz.+DT:/DS:ADSP\\\tDT:/g" -)

echo "`date`: Realigning sample $sampleName" 
TMP_DIR="$sampleName"
mkdir -p $TMP_DIR

# Will use .bam for intermediate files regardless of the specified output
# format (for simplicity).
sample_sorted_by_name="${TMP_DIR}/${sampleName}.name_sorted.bam"

echo "`date`: Sample data sorted by name will be written to $sample_sorted_by_name"

# Sorting the sample file by name to extract read pairs for realignment
echo "`date`: Sorting $sample_input by name"
if ! time samtools sort -@ $n_threads -n -m ${mem_per_thread}G $sample_input > $sample_sorted_by_name; then
	echo "sort failed"
	exit 1
fi

# Defining the fastq names
tmp_fq1="${TMP_DIR}/tmp_${sampleName}_R1.fastq"
tmp_fq2="${TMP_DIR}/tmp_${sampleName}_R2.fastq"


# Export CRAM_REFERENCE so bamtofastq has it when the input is a .cram file
export CRAM_REFERENCE=$original_ref

# Extracting all reads from the bam/cram and converting to fastq for
# realignment
echo "`date`: Converting to Fastq:"
if ! time bedtools bamtofastq -i $sample_sorted_by_name -fq $tmp_fq1 -fq2 $tmp_fq2; then
	echo "bamtofastq failed"
	exit 1
fi

# Defining the fastq names
fq1="${TMP_DIR}/${sampleName}_R1.fastq"
fq2="${TMP_DIR}/${sampleName}_R2.fastq"

# Fix bad formatting this sample's read names (uses old format that makes bwa choke). This was used
# for one specific sample from dbGap that was originally aligned many years ago and had some funky
# read names that bwa choked on. I'm leaving this here just in case it's ever useful. This removes
# the '/1' and '/2' entirely.
#
# The variable names above and below ($fq1 and $fq2) would need to be modified to plug this back in
# to the script. Otherwise will overwrite $fq1 and $fq2 with nothing.
# awk '{ if(NR%4==1) print gensub(/(^@[[:digit:]]*)(\/[[:digit:]]$)/, "\\1", "g", $0); else print $1}' $tmp_fq1 > $fq1
# awk '{ if(NR%4==1) print gensub(/(^@[[:digit:]]*)(\/[[:digit:]]$)/, "\\1", "g", $0); else print $1}' $tmp_fq2 > $fq2

tmp_sample="${TMP_DIR}/${sampleName}.unsorted.${out_format}"

# Aligning fastqs to $align_to_ref and converting to .bam file

if [ "$out_format" = "bam" ]; then
	echo "`date` Aligning to $sampleName to ${align_to_ref}. Output will be .bam file"
	time bwa mem -M -R "$RG" -t $n_threads $align_to_ref $fq1 $fq2 | samtools view -hb - > $tmp_sample
else
	echo "`date` Aligning to $sampleName to ${align_to_ref}. Output will be .cram file"
	time bwa mem -M -R "$RG" -t $n_threads $align_to_ref $fq1 $fq2 | samtools view -C -T $align_to_ref - > $tmp_sample
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


final_sample_output_file="${sampleName}.${ref_tag}.${out_format}"

# Sorting the aligned sample file
echo "`date` Sorting final sample $out_format"
if ! time samtools sort -@ $n_threads -m ${mem_per_thread}G $tmp_sample > $final_sample_output_file; then
	echo "ERROR: Sorting final sample $out_format failed"
	exit 1
fi

# Indexing the sorted .bam/cram
echo "`date` Indexing final $out_format file ($final_sample_output_file)"
if ! time samtools index $final_sample_output_file; then
	echo "ERROR: Indexing final $out_format failed"
	exit 1
fi


# Validating the final .bam/cram for output.
echo "`date` Validating $out_format"
if ! time samtools quickcheck $final_sample_output_file; then
	echo "quickcheck failed"
	exit 1
else
	echo "quickcheck happy!"
fi


echo "`date` DONE!"
#rm -rfv $TMP_DIR
