process CREATE_BED_FILE {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/CREATE_BED_FILE", mode: 'copy')
	// publishDir "${projectDir}/../results/create_bed_file", mode: 'copy'

	label 'CREATE_BED_FILE'

	input:
		path(low_depth)
		path(low_mapq)
		file(ref)
		file(ref_index)
		path(annotation_bed)
		val(sequencer)
		val(ref_tag)
		// val(threads)

	output:
		// path 'illuminaRL100.hg38*' 
		path "${sequencer}.${ref_tag}*" 

	script:
	"""
	bash camo_gene_pipeline.sh \
		-d ${low_depth} \
		-m ${low_mapq} \
		-g ${ref} \
		-a ${annotation_bed} \
		-s ${sequencer} \
		-v ${ref_tag} \
		-t $task.cpus \
	"""
}
