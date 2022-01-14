
workflow RESCUE_CAMO_VARS_WF {
    take:
        realign_bed
        gatk_bed
        mask_ref_prefix
        bam_path

    main:
        PREP_BATCH_FILES(bam_path)

        // bam_sets_ch = Channel.fromPath("$PREP_BATCH_FILES.out.bam_set/bam_set_*", checkIfExists: true)
        // RESCUE_CAMO_VARS_PROC(realign_bed, gatk_bed, mask_ref_prefix, bam_sets_ch)
        RESCUE_CAMO_VARS_PROC(realign_bed, gatk_bed, mask_ref_prefix, PREP_BATCH_FILES.out.bam_sets.flatten())
}


/*
 * The process to rescue camouflaged variants
 */
process RESCUE_CAMO_VARS_PROC {
	
	label 'RESCUE_CAMO_VARS_PROC'

	input:
		val(realign_bed)
		val(gatk_bed)
		val(mask_ref_prefix)
        val bam_set

// 	output:
// 		path '*.g.vcf', emit: camo_gvcfs

	script:
	"""
		bash rescue_camo_vars.sh $realign_bed $gatk_bed $mask_ref_prefix $bam_set
	"""
}

/*
 * The process to prepare batch files for the rescue process. This
 * will split the set of samples into file lists of 50 sam/bam/cram
 * files.
 */
process PREP_BATCH_FILES {

	label 'local'

	input:
		val(bam_path)

	output:
		path 'bam_set_*', emit: bam_sets

	script:
	"""
    # The loop below will produce errors if files with any of the specified
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
