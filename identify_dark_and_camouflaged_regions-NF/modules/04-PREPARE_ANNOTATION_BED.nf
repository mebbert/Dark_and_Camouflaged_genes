process PREPARE_ANNOTATION_BED_PROC {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-PREPARE_ANNOTATION_BED", mode: 'copy')

	label 'local'

	input:
		path(annotation_bed) //path to gene annotation gff3

	output:
		path '*.annotation.bed', emit: prepped_anno_bed
		path 'intermediate_annotation_files', emit: intermediate_files

	script:
		"""
		bash create_annotation_bed.sh ${annotation_bed}
		"""
}
