#!/bin/usr bash

# Enable bash debugging to log all commands
set -x

echo "Job running on: `hostname`"

# The masked ref fasta
ref=$1

# the current ref fasta to be used for bamtofastq for crams
current_ref=$2

# The bed file to use for extracting reads from each sample
extraction_bed=$3		
echo "Extraction bed: $extraction_bed"

# The bed file to use for calling variants. This will make GATK much faster.
# Using a different bed file here because the regions are restricted to the
# exact camo CDS regions.
gatk_bed=$4
echo "GATK bed: $gatk_bed"


# A file with the list of bams to run
bam_list=$5
echo "Bam list file: $bam_list"

# The max number of repeat regions to consider
max_repeat=$6
echo "Max repeats: $max_repeat"

# The min number of repeat regions to consider. By definition, two is the fewest.
# Otherwise it wouldn't be camouflaged.
min_repeat=2 # 

# Whether to clean tmp files on the fly. This process could could overwhelm a
# file system (depending on how large the run is).
clean_tmp_files=$7

# Create a new list for GATK for each ploidy
for repeat_num in $(seq $min_repeat 1 $max_repeat)
do
	ploidy=$(( 2 * $repeat_num ))
	ploidy_list="gvcfs_ploidy${ploidy}.list"
	gvcf_list[$repeat_num]=$ploidy_list
done 

# Regex to extract sample name from Bam header
sample_regex="SM:([A-Za-z0-9_\-]+)"

threads=$8

hardcoded_ploidy=$9

# Create GATK environment variable for convenience.
GATK="gatk --java-options -Xmx20G "

# Track first and last sample to provide unique name to final output .gvcf for this
# batch. Cannot use JobID because we don't know which scheduler (if any) will be
# used. Cannot use seconds since Epoch because two jobs might finish at the same time.
index=0
first_sample=""
last_sample=""

# For each input bam, extract reads from camouflaged regions, align
# to masked genome, and call variants with GATK.
while read bam
do

	filename=$(basename -- "$bam")
	extension="${filename##*.}"
	
	echo "Extension: $extension"

	# Determine if we're working with .bam or .cram (',,' == convert $extension to all lower case)
	# for "case-insensitive" comparison.
	is_bam=false
	if [[ ${extension,,} == "bam" ]]; then
		is_bam=true
	elif [[ ${extension,,} == "cram" ]]; then
		is_bam=false;
	else
		echo "[rescue_camo_variants.sh]: Please convert the following file to either a .bam or .cram: $bam"
		exit 1
	fi


	# Extract the sample name from .bam's read group header
	RG=$(lsamtools view -H $bam | grep '^@RG' | tail -1)
	[[ $RG =~ $sample_regex ]]
	sampleName=${BASH_REMATCH[1]}

	echo "Found sample name: $sampleName"
	if [[ "$index" -eq 0 ]]; then
		first_sample=$sampleName
	fi

	# Increment
	index=$((index+1))

	# Always set last_sample to sampleName. It will be on the last loop.
	last_sample=$sampleName

	tmp_bam=${sampleName}.sorted_by_name.bam


	#######################################
	# Rescue variants for each ploidy set #
	#######################################
	for repeat_num in $(seq $min_repeat 1 $max_repeat)
	do
		# Define the ploidy as 2 * the number of identical regions
		ploidy=$(( 2 * $repeat_num ))

		# Concat all .bed regions for a given repeat set into a single "regions" list that can be passed to
		# samtools on the command line (not as .bed file). For some reason, samtools does not use the .bam index
		# if you pass a .bed file using the -L <bed> option (at least at the time
		# this was written. 
		#
		# As an example, when \$repeat_num == 2, this command will only provide the .bed entries
		# that have only two identical regions.
		#
		# The fist part is skipping comment lines. 'paste' is a unix command that pastes file lines together.
		regions=$(awk "\$5 == $repeat_num { print \$1\":\"\$2\"-\"\$3 }" $extraction_bed | paste -sd " ")
		if [[ -z $regions ]]; then continue; fi #if regions is empty skip 

		# Extract only reads with MapQ < 10
		time lsamtools view -h $bam $regions | \
			awk '$5 < 10 || $1 ~ "^@"' | \
			lsamtools view -hb - | \
			lsamtools sort -n -m 16G -o $tmp_bam -


		####################
		# Convert to fastq #
		####################

		# Define tmp fastq files
		fq1=${sampleName}_R1.fastq
		fq2=${sampleName}_R2.fastq

		# TODO: still need to check if these work and will fix our problem
		#if [[ "$(head -1 $fq1)" =~ ^@\\1.* ]]; then
		#	sed -i 1,4d $fq1
		#fi
		#if [[ "$(head -1 $fq2)" =~ ^@\\2.* ]]; then
		#	sed -i 1,4d $fq2
		#fi

		# TODO: Handle bam & cram!
		export CRAM_REFERENCE=$current_ref
		time bedtools bamtofastq -i $tmp_bam -fq $fq1 -fq2 $fq2 2> /dev/null

		##################
		# Align with BWA #
		##################
		aligned_sam=${sampleName}_ploidy_${ploidy}.sam
		final_bam=${sampleName}.ploidy_${ploidy}.bam
		# RG="@RG\tID:group1\tSM:$sampleName\tPL:illumina\tLB:lib1\tPU:unit1"

		# Replace literal tab characters with '\t' so bwa doesn't get angry. The
		# double '//' means global search/replace.
		RG=${RG//$'\t'/"\t"}
		time bwa mem -M  \
			-R "$RG" \
			-t $threads \
			$ref \
			$fq1 $fq2 > $aligned_sam

		# Sort and index bam
		time lsamtools view -bt $ref $aligned_sam | lsamtools sort -@ $threads -m 6G -o $final_bam -
		lsamtools index $final_bam


		##################
		# Call mutations #
		##################

		# I don't see how we can perform base recalibrator since
		# these regions have never been characterized.
		gvcf=${sampleName}.ploidy_${ploidy}.g.vcf
		echo $gvcf >> ${gvcf_list[$repeat_num]}
		
		# For clarity, we previously explicitly invoked '--genotyping-mode DISCOVERY' (GATK v3.X)
		# even though it was the default setting. That setting has been removed in GATK v4.X but
		# apparently remains the default: 
		# https://gatk.broadinstitute.org/hc/en-us/community/posts/360077648352-GATK-4-2-0-0-Haplotype-caller-genotyping-mode-DISCOVERY-Option-
		time $GATK HaplotypeCaller \
			--native-pair-hmm-threads "${threads}" \
			-R $ref \
			-I $final_bam \
			-L $gatk_bed \
			--sample-ploidy "${hardcoded_ploidy}" \
			--emit-ref-confidence GVCF \
			--dont-use-soft-clipped-bases \
			-O $gvcf

		# Cleaning up because this pipeline may create thousands of files, depending 
		# on how large the data set is.
		if [ "$clean_tmp_files" = true ] ; then
			rm ${aligned_sam}* ${tmp_bam}* ${final_bam}* $fq1 $fq2
		fi
	done
done < $bam_list


#####################################
# Combine gVCFs for each ploidy set #
#####################################
for repeat_num in $(seq $min_repeat 1 $max_repeat)
do

	# If this gvcf_list does not exist, skip
	if [ ! -f "${gvcf_list[$repeat_num]}" ]; then
		continue
	fi

	ploidy=$(( 2 * $repeat_num ))
	comb_gvcf_file=camo_batch_gvcfs/ploidy_${ploidy}/samples_${first_sample}_through_${last_sample}.g.vcf

	# Create a results directory for this ploidy since we know
	# there will be results for it.
	mkdir -p camo_batch_gvcfs/ploidy_${ploidy}/

	time $GATK CombineGVCFs \
		-R $ref \
		-O $comb_gvcf_file \
		--variant ${gvcf_list[$repeat_num]}

	# Index the gvcf so that it can be accessed 'randomly'
	# by CombineGVCFs
	time $GATK IndexFeatureFile \
		-I $comb_gvcf_file

done

# Cleaning up because this pipeline may create thousands of files, depending 
# on how large the data set is.
if [ "$clean_tmp_files" = true ] ; then
	rm *.vcf*
fi
echo "COMMAND(`date`): DONE."
