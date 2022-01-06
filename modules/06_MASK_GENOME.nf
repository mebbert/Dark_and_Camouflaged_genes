process step_06_MASK_GENOME {
	
	label 'local'

	input:
		path(align_to_bed)
		path(ref)
		path(ref_idx)
		val(mask_ref_prefix)

	output:
		path '*.fa', emit: masked_fastas
		path '*.dict', emit: sequence_dicts
		path '*.fa.*', emit: fasta_idx 

	script:
		"""
		bash create_camo_mask.sh $align_to_bed $ref $mask_ref_prefix

		"""
}
