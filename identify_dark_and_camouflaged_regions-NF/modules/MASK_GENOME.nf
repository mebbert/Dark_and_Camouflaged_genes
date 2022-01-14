process MASK_GENOME {

    publishDir 'results/MASK_GENOME'
	
	label 'MASK_GENOME'

	input:
		path(align_to_bed)
		path(ref)
		path(ref_idx)
		val(mask_ref_prefix)

	output:
		path('*.fa', emit: masked_ref_fasta)
		path('*.fa.*', emit: masked_ref_idx)

	script:
		"""
		bash create_camo_mask.sh $align_to_bed $ref $mask_ref_prefix
		"""
}
