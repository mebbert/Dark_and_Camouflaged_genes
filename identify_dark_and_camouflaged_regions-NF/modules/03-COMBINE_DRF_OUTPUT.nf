
/*
 * Combine DRF output from all samples, and then separate the results into
 * separate 'dark-by-depth' and 'dark-by-MAPQ' .bed files. The
 * 'result_prefix' string is used for naming final results output.
 */
process COMBINE_DRF_OUTPUT_PROC {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/03-COMBINE_DRF_OUTPUT", mode: 'copy')

	label 'COMBINE_DRF_OUTPUT_PROC'

	input:
		path(low_mapq_bed_list)
		val(result_prefix)

	output:
		path('*.dark.low_depth.bed.gz'), emit: low_depth_out
		path('*.dark.low_mapq.bed.gz'), emit: low_mapq_out

	script:
	"""
	bash combine_DRF.sh '${low_mapq_bed_list}' $result_prefix $task.cpus
	"""
}
