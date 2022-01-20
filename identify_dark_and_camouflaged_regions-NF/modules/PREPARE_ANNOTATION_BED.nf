process PREPARE_ANNOTATION_BED {
	
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
