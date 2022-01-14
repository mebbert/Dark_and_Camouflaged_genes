process step_05_CREATE_BED_FILE {
	
	publishDir '/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/results/step_05', mode: 'copy'

	label 'step_05'

	input:
		path(low_depth)
		path(low_mapq)
		file(ref)
		file(ref_index)
		path(annotation_bed)
		val(sequencer)
		val(ref_tag)
		val(threads)
		val(output_dir)

	output:
		path 'illuminaRL100.hg38*' 

	script:
	"""
	bash camo_gene_pipeline.sh \
		-d ${low_depth} \
		-m ${low_mapq} \
		-g ${ref} \
		-a ${annotation_bed} \
		-s ${sequencer} \
		-v ${ref_tag} \
		-t ${threads} \
		-r ${output_dir}
	"""
}
