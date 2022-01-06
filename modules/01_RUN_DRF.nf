process step_01_RUN_DRF {
	
	label 'step_01'

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
