process MASK_GENOME {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/MASK_GENOME", mode: 'copy')
	
	label 'MASK_GENOME'

	input:
		path(align_to_bed)
		path(ref)
		path(ref_idx)
		val(mask_ref_prefix)

	output:
		path('*.fa', emit: masked_ref_fasta)
		path('*.fa.*', emit: masked_ref_idx)
		path('*.dict', emit: masked_ref_idx)

	script:
		"""
		bash create_camo_mask.sh $align_to_bed $ref $mask_ref_prefix
		"""
}
