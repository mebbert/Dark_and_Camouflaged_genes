process MASK_GENOME_PROC {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/06-MASK_GENOME", mode: 'copy')
	
	label 'MASK_GENOME_PROC'

	input:
		tuple val(sample_name), path(align_to_bed)
		val(mask_ref_prefix)

	output:
		path('*.fa', emit: masked_ref_fasta)
		path('*.fa.*', emit: masked_ref_idx)
		path('*.dict', emit: masked_ref_dict)

	script:
		"""
		bash create_camo_mask.sh $align_to_bed $params.align_to_ref $mask_ref_prefix
		"""
}
