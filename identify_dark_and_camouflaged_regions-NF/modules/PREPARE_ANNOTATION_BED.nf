process PREPARE_ANNOTATION_BED {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/PREPARE_ANNOTATION_BED", mode: 'copy')

	label 'local'

	input:
		path(annotation_bed) //path to gene annotation gff3

	output:
		path '*.bed', emit: prepped_anno_bed

	script:
		"""
		bash create_annotation_bed.sh ${annotation_bed}
		"""
}
