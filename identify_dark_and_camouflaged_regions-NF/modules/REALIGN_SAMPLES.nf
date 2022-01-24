process REALIGN_SAMPLES {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/REALIGN_SAMPLES", mode: 'copy')

	label 'REALIGN_SAMPLES'

	input:
		path(cram)
		val(align_to_ref)
		val(ref_tag)
		val(original_ref)
		
	output:
		path '*.bam', emit: final_bams

	script:
	"""
	bash realign_bwa.sh ${cram} ${align_to_ref} ${ref_tag} ${original_ref}
	"""
}
