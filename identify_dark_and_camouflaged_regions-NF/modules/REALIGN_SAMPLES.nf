process REALIGN_SAMPLES {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/REALIGN_SAMPLES", mode: 'copy')

	label 'REALIGN_SAMPLES'

	input:
		path(sample_input_file)
		val(align_to_ref)
		val(ref_tag)
		val(original_ref)
        val(output_format)


	output:
        path '*.{bam,cram}', emit: final_alignment
        path '*.{bam,cram}.*', emit: final_alignment_index

	script:

    def avail_mem = task.memory ? task.memory.toGiga() : 0

    /*
     * Calculate the mem per thread. Divide by an extra thread to provide
     * some buffer.
     */
    def mem_per_thread = (avail_mem).intdiv(task.cpus + 1)

	"""
    echo "Sample file: $sample_input_file"
	bash realign_bwa.sh ${sample_input_file} ${align_to_ref} ${ref_tag} ${original_ref} ${output_format} $task.cpus $mem_per_thread
	"""
}
