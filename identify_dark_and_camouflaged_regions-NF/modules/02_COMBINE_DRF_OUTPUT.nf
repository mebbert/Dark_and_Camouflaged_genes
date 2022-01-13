process step_02_COMBINE_DRF_OUTPUT {
	
	label 'step_02'

	input:
		file(low_mapq_bed)
		val(result_prefix)

	output:
		path '*.dark.low_depth.bed', emit: low_depth_out
		path '*.dark.low_mapq.bed', emit: low_mapq_out

	script:
	"""
	echo $low_mapq_bed
	echo $result_prefix
	bash combine_DRF.sh '${low_mapq_bed}' $result_prefix
	"""
}
