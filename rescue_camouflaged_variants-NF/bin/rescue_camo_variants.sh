#!/bin/usr bash

# Enable bash debugging to log all commands
set -x

echo "Job running on: `hostname`"

# The masked ref fasta
ref=$1

# The bed file to use for extracting reads from each sample
extraction_bed=$2		
echo "Extraction bed: $extraction_bed"

# The bed file to use for calling variants. This will make GATK much faster.
# Using a different bed file here because the regions are expanded a bit
# since reads extend past.
gatk_bed=$3
echo "GATK bed: $gatk_bed"

echo "###################################"
echo "PATH: $PATH"
echo "###################################"

# Create GATK environment variable for convenience.
GATK="gatk --java-options -Xmx20G"

# A file with the list of bams to run
bam_list=$4
echo "Bam list file: $bam_list"

# The max number of repeat regions to consider
max_repeat=$5
echo "Max repeats: $max_repeat"

# The min number of repeat regions to consider. By definition, two is the fewest.
# Otherwise it wouldn't be camouflaged.
min_repeat=2 # 


# Create a new list for GATK for each ploidy
for repeat_num in $(seq $min_repeat 1 $max_repeat)
do
	ploidy=$(( 2 * $repeat_num ))
	mkdir -p camo_gvcfs/ploidy_${ploidy}/
	ploidy_list="gvcfs_ploidy${ploidy}.list"
	gvcf_list[$repeat_num]=$ploidy_list
done 

# Regex to extract sample name from Bam header
sample_regex="SM:([A-Za-z0-9_\-]+)"

threads=1

# Track first and last sample to provide unique name to final output .gvcf. Cannot
# use JobID because we don't know which scheduler (if any) will be used. Cannot use
# seconds since Epoch because two jobs might finish at the same time.
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
	RG=$(samtools view -H $bam | grep '^@RG' | tail -1)
	[[ $RG =~ $sample_regex ]]
	sampleName=${BASH_REMATCH[1]}

	echo "Found sample name: $sampleName"
	if [[ "$index" -eq 0 ]]; then
		first_sample=$sampleName

	# Always set last_sample to sampleName. It will be on the last loop.
	last_sample=$sampleName

	tmp_bam=${sampleName}.sorted.bam

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
		time samtools view -h $bam $regions | \
			awk '$5 < 10 || $1 ~ "^@"' | \
			samtools view -hb - | \
			samtools sort -n -m 16G -o $tmp_bam -


		####################
		# Convert to fastq #
		####################

		# Define tmp fastq files
		fq1=${sampleName}_R1.fastq
		fq2=${sampleName}_R2.fastq

		# TODO: Handle bam & cram!
		time bedtools bamtofastq -i $tmp_bam -fq $fq1 -fq2 $fq2 2> /dev/null

		##################
		# Align with BWA #
		##################
		aligned_sam=${sampleName}.sam
		final_bam=${sampleName}.ploidy_${ploidy}.bam
		# RG="@RG\tID:group1\tSM:$sampleName\tPL:illumina\tLB:lib1\tPU:unit1"

		# Replace literal tab characters with \t so bwa doesn't get angry. The
		# double '//' means global search/replace.
		RG=${RG//$'\t'/"\t"}
		time bwa mem -M  \
			-R "$RG" \
			-t $threads \
			$ref \
			$fq1 $fq2 > $aligned_sam

		# Sort and index bam
		time samtools view -bt $ref $aligned_sam | samtools sort -@ $threads -m 6G -o $final_bam -
		samtools index $final_bam


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
			-R $ref \
			-I $final_bam \
			-L $gatk_bed \
			--sample-ploidy $ploidy \
			--emit-ref-confidence GVCF \
			--dont-use-soft-clipped-bases \
			-O $gvcf

		# Cleaning up because this pipeline may create thousands of files, depending 
		# on how large the data set is.
		# rm ${aligned_sam}* ${tmp_bam}* ${final_bam}* $fq1 $fq2
	done
done < $bam_list

#####################################
# Combine gVCFs for each ploidy set #
#####################################
for repeat_num in $(seq $min_repeat 1 $max_repeat)
do
	ploidy=$(( 2 * $repeat_num ))
	comb_gvcf_file=camo_gvcfs/ploidy_${ploidy}/${first_sample}_to_${last_sample}.g.vcf
	time $GATK CombineGVCFs \
		-R $ref \
		-O $comb_gvcf_file \
		--variant ${gvcf_list[$repeat_num]}

done

# Cleaning up because this pipeline may create thousands of files, depending 
# on how large the data set is.
# rm *.vcf*
echo "COMMAND(`date`): DONE."
