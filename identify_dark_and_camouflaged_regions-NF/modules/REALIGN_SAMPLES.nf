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
        val(output_format)


	output:
        path '*.{bam,cram}*', emit: final_alignments

	script:

    def avail_mem = task.memory ? task.memory.toGiga() : 0

    /*
     * Calculate the mem per thread. Divide by an extra thread to provide
     * some buffer.
     */
    def mem_per_thread = (avail_mem).intdiv(task.cpus + 1)

	"""
	bash realign_bwa.sh ${cram} ${align_to_ref} ${ref_tag} ${original_ref} ${output_format} $task.cpus $mem_per_thread
	"""
}
