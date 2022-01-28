process RUN_DRF {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/RUN_DRF", mode: 'copy')

	label 'RUN_DRF'

	input:
		path(sample_input_file)
		val(ref)
		path(jar)

	output:
		path '*.dark.low_mapq.bed', emit: low_mapq_bed

	script:
	"""
	bash run_DRF.sh ${sample_input_file} ${ref} ${jar}
	"""
}
