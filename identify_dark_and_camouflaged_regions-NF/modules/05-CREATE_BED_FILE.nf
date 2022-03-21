process CREATE_BED_FILE_PROC {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/05-CREATE_BED_FILE", mode: 'copy')

	label 'CREATE_BED_FILE_PROC'

	input:
		path(low_depth)
		path(low_mapq)
		path(annotation_bed)

	output:
		path("*mask_bed.sorted.bed"), emit: mask_bed
		path "*" 

	script:
	"""
	bash camo_gene_pipeline.sh \
		-d ${low_depth} \
		-m ${low_mapq} \
		-g ${params.align_to_ref} \
		-a ${annotation_bed} \
		-s ${params.sequencer_tag} \
		-v ${params.align_to_ref_tag} \
		-t $task.cpus \
	"""
}
