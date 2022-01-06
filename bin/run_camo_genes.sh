#!/bin/usr bash

TEMPORARY_DIR="tmp/$JOB_ID"

#	set up function.  this isn't called/run here.  It's only used
#   if the job is canceled via a signal
cleanup_scratch()
{
	echo "Deleting inside signal handler, meaning I probably either hit the walltime, or deleted the job"

	#remove the remaining data in $TEMPORARY_DIR
	#rm -rfv "$TEMPORARY_DIR"
	echo "---"
	echo "Signal handler ending time:"
	date
	exit 0
}


# Associate the function "cleanup_scratch" with the TERM signal, which is usually how jobs get killed
trap 'cleanup_scratch' TERM


echo "Job running on: `hostname`"


echo "Creating TMP dir: $TEMPORARY_DIR"
mkdir -p $TEMPORARY_DIR

# The bed file to use for extracting reads from each sample
realign_bed=$1		
echo $realign_bed

# The bed file to use for calling variants. This will make GATK much faster.
# Using a different bed file here because the regions are expanded a bit
# since reads extend past.
gatk_bed=$2			

threads=1

#Path prefix to camo-masked references
# Should look like: path/to/refs/b37-camo_mask
ref_prefix=$3
filtered_bam_list=$4
DATA_DEST_DIR='/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/final_files'

# Regex to extract sample name from Bam header
regex="SM:([A-Za-z0-9_\-]+)"

max_repeat=5
for repeat_num in $(seq 2 1 $max_repeat)
do
	ploidy=$(( $repeat_num * 2 ))
	mkdir -p ${DATA_DEST_DIR}/camo_gvcfs/ploidy_${ploidy}/
	ploidy_list="$TEMPORARY_DIR/gvcfs_ploidy${ploidy}.list"
	gvcf_list[$repeat_num]=$ploidy_list
done

while read bam
do

	RG=$(samtools view -H $bam | grep '^@RG' | tail -1)
	[[ $RG =~ $regex ]]
	sampleName=${BASH_REMATCH[1]}

	echo "Sample name: $sampleName"

	tmp_bam=$TEMPORARY_DIR/${sampleName}.sorted.bam

	for repeat_num in $(seq 2 1 $max_repeat)
	do
		ploidy=$(( 2 * $repeat_num ))
		echo "testing ploidy level $ploidy"
		ref="${ref_prefix}.ploidy_${ploidy}.fa"

		# concat bed regions into a single "regions" list. Samtools does not use the .bam index
		# if you use the -L <bed> option.
		#
		# The fist part is skipping comment lines. 'paste' is a unicx command that pastes file lines together.
		regions=$(awk "\$5 == $repeat_num { print \$1\":\"\$2\"-\"\$3 }" $realign_bed | paste -sd " ")
		if [[ -z $regions ]]; then continue; fi #if regions is empty skip 

		# Extract reads from regions and sort by read name.
		#echo "samtools view -hb $bam | samtools sort -n -m 10G -o $tmp_bam"
		#samtools view -hb $bam $regions | samtools sort -@ $threads -n -m 10G -o $tmp_bam

		#Extract only reads with MapQ < 10
		echo "samtools view -h $bam regions | awk '\$5 < 10' | samtools view -hb > $tmp_bam"
		time samtools view -h $bam $regions | \
			awk '$5 < 10 || $1 ~ "^@"' | \
			samtools view -hb - | \
			samtools sort -n -m 16G -o $tmp_bam -

		# Convert to fastq
		fq1=$TEMPORARY_DIR/${sampleName}_R1.fastq
		fq2=$TEMPORARY_DIR/${sampleName}_R2.fastq

		echo "bedtools bamtofastq -i $tmp_bam -fq $fq1 -fq2 $fq2"
		time bedtools bamtofastq -i $tmp_bam -fq $fq1 -fq2 $fq2 2> /dev/null

		# Remove tmp_bam
		rm ${tmp_bam}*

		# Align with BWA
		aligned_sam=$TEMPORARY_DIR/${sampleName}.sam
		final_bam=$TEMPORARY_DIR/${sampleName}.ploidy_${ploidy}.bam
		RG="@RG\tID:group1\tSM:$sampleName\tPL:illumina\tLB:lib1\tPU:unit1"
		echo "bwa mem -M -R $RG -t $threads $ref $fq1 $fq2 > ${aligned_sam}"
		time bwa mem -M  \
			-R $RG \
			-t $threads \
			$ref \
			$fq1 $fq2 > $aligned_sam

		# Convert to sorted bam
		echo "samtools view -bt $ref $aligned_sam | samtools sort -@ $threads -m 6G -o $final_bam -"
		time samtools view -bt $ref $aligned_sam | samtools sort -@ $threads -m 6G -o $final_bam -

		# Remove sam and fastqs
		rm $aligned_sam
		rm $fq1
		rm $fq2

		# index bam
		echo "samtools index $final_bam"
		samtools index $final_bam


		# Call mutations
		# I don't see how we can perform base recalibrator since
		# these regions have never been characterized.
		gvcf=$TEMPORARY_DIR/${sampleName}.${repeat_num}.g.vcf

		echo $gvcf >> ${gvcf_list[$repeat_num]}
		
		time gatk HaplotypeCaller \
			-R $ref \
			-I $final_bam \
			-L $gatk_bed \
			--sample-ploidy $ploidy \
			--genotyping-mode DISCOVERY \
			--emit-ref-confidence GVCF \
			--dont-use-soft-clipped-bases \
			-O $gvcf

		rm ${final_bam}*
	done
done < $filtered_bam_list

for repeat_num in $(seq 2 1 $max_repeat)
do
	ploidy=$(( $repeat_num * 2 ))
	ref="${ref_prefix}.ploidy_${ploidy}.fa"
	comb_gvcf_file=$DATA_DEST_DIR/camo_gvcfs/ploidy_${ploidy}/$JOB_ID.g.vcf
	echo "$GATK -T CombineGVCFs -R $ref -o $comb_gvcf_file  --variant ${gvcf_list[$repeat_num]}"
	gatk CombineGVCFs \
		-R $ref \
		-O $comb_gvcf_file \
		--variant ${gvcf_list[$repeat_num]}
done

#############
## CLEAN UP #
#############
#change to a safe location
#cd "$HOME"

#remove the data in $TEMPORARY_DIR
echo "Removing tmp files files..."
rm -rfv "$TEMPORARY_DIR"
echo "---"

echo "COMMAND(`date`): DONE."
