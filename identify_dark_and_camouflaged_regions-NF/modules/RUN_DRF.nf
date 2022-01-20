process RUN_DRF {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/RUN_DRF", mode: 'copy')

	label 'run_drf'

	input:
		path(bam)
		file(ref)
		file(ref_index)
		path(jar)

	output:
		path '*.dark.low_mapq.bed', emit: low_mapq_bed

	script:
	"""
	bash run_DRF.sh ${bam} ${ref} ${jar}
	"""
}
