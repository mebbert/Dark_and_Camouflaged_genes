process REALIGN_BAMS {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/REALIGN_BAMS", mode: 'copy')

	label 'realign_bams'

	input:
		path(cram)
		file(ref)
		val(ref_tag)
		file(cram_ref)
		file(ref_ind)
		file(ref_ind_amb)
		file(ref_ind_pac)
		file(ref_ind_ann)
		file(ref_ind_sa)
		file(ref_ind_bwt)

	output:
		path '*.bam', emit: final_bams

	script:
	"""
	bash realign_bwa.sh ${cram} ${ref} ${ref_tag} ${cram_ref}
	"""
}
