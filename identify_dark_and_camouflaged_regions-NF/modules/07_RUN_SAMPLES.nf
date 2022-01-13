process REALIGN_SAMPLES_PROC {
	
	label 'REALIGN_SAMPLES_PROC'

	input:
		val(realign_bed)
		val(gatk_bed)
		val(mask_ref_prefix)
		val(bam_set_ch)

	output:
		// path '*.g.vcf', emit: camo_gvcfs

	script:
	"""
        echo "########################"
		echo "Filtered list: $bam_set"
        echo "########################"
        echo ""
		bash run_camo_genes.sh $realign_bed $gatk_bed $mask_ref_prefix $bam_set
	"""
}

process step_07_PREP {

	label 'local'

	input:
		val(bam_path)

	output:
		path 'bam_set*', emit: bam_set

	script:
	"""

    # The loop below will produce errors if no files with any of the specified
    # extensions aren't found (e.g., if there are no *.sam files). Set the
    # `nullglob` to hide errors. Also set `nocaseglob` to ignore case
    # for the file extension (.e.g., *.BAM, *.BaM, etc.).
    shopt -s nullglob # Sets nullglob
    shopt -s nocaseglob # Sets nocaseglob

    # Loop over all files ending with .sam, .bam, and .cram and print
    # to file. Then split file into sets of 50 so we can process sets of
    # 50 .bam files per job submission.

	for bam in $bam_path/*.{sam,bam,cram}
	do
		echo \$bam
	done > bam_list.txt

    # Split bam_list.txt into sets of 50 bam files with prefix 'bam_set_'
    # using numeric suffixes (-d; e.g., 00, 01, etc.).
	split -dl 50 bam_list.txt 'bam_set_'

    # Unset nullglob and nocaseglob
    shopt -u nullglob # Unsets nullglob
    shopt -u nocaseglob # Unsets nocaseglob
	"""
}
