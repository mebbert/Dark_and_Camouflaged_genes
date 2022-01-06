process step_07_RUN_SAMPLES {
	
	label 'step_07'

	input:
		val(realign_bed)
		val(gatk_bed)
		val(mask_ref_prefix)
		path(filtered_list)

	output:
		path '*.g.vcf', emit: camo_gvcfs

	script:
	"""
		echo $filtered_list
		bash run_camo_genes.sh $realign_bed $gatk_bed $mask_ref_prefix $filtered_list
	"""
}

process step_07_PREP {

	label 'local'

	input:
		val(filtered_bam_list)

	output:
		path 'split_list*', emit: filtered_lists

	script:
	"""

	for file_bam in $filtered_bam_list
	do
		echo \$file_bam
	done > files_to_split.txt

	split -l 50 files_to_split.txt 'split_list_'
	"""
}
