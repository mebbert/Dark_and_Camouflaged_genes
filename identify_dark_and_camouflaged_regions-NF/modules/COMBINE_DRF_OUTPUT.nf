process COMBINE_DRF_OUTPUT {
	
	label 'combine_drf_output'

	input:
		file(low_mapq_bed_list)
		val(result_prefix)

	output:
		path '*.dark.low_depth.bed', emit: low_depth_out
		path '*.dark.low_mapq.bed', emit: low_mapq_out

	script:
	"""
	echo $low_mapq_bed_list
	echo $result_prefix
	bash combine_DRF.sh '${low_mapq_bed_list}' $result_prefix
	"""
}
