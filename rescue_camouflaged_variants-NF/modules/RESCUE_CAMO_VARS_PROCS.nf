import groovy.io.FileType

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

workflow RESCUE_CAMO_VARS_WF {
    take:
        bam_path
        extraction_bed
        masked_ref_fastas
        gatk_bed
        n_samples_per_batch
        max_repeats_to_rescue

    main:

        println bam_path
        PREP_BATCH_FILES_PROC(bam_path, n_samples_per_batch)
        RESCUE_CAMO_VARS_PROC(
                                extraction_bed,
                                masked_ref_fastas,
                                gatk_bed,
                                PREP_BATCH_FILES_PROC.out.bam_sets.flatten(),
                                max_repeats_to_rescue
                              )

        GENERATE_MASTER_GVCF_LIST_PROC(RESCUE_CAMO_VARS_PROC.camo_gvcf_dirs.collect())
}


/*
 * The process to prepare batch files for the rescue process. This
 * will split the set of samples into file lists of 50 sam/bam/cram
 * files.
 */
process PREP_BATCH_FILES_PROC {

	label 'local'

	input:
		val(bam_path)
        val(n_samples_per_batch)

	output:
		path('bam_set_*', emit: bam_sets)

	script:
	"""
    # The loop below will produce errors if files with any of the specified
    # extensions aren't found (e.g., if there are no *.sam files). Set the
    # `nullglob` to hide errors. Also set `nocaseglob` to ignore case
    # for the file extension (.e.g., *.BAM, *.BaM, etc.).
    shopt -s nullglob # Sets nullglob
    shopt -s nocaseglob # Sets nocaseglob

    # Loop over all files ending with .sam, .bam, and .cram and print
    # to file. Then split file into sets of $n_samples_per_batch so we can
    # process sets of $n_samples_per_batch .bam files per job submission.

	for bam in $bam_path/*.{bam,cram}
	do
		echo \$bam
	done > bam_list.txt

    # Split bam_list.txt into sets of 50 bam files with prefix 'bam_set_'
    # using numeric suffixes (-d; e.g., 00, 01, etc.).
	split -dl $n_samples_per_batch bam_list.txt 'bam_set_'

    # Unset nullglob and nocaseglob
    shopt -u nullglob # Unsets nullglob
    shopt -u nocaseglob # Unsets nocaseglob
	"""
}


/*
 * The process to rescue camouflaged variants. Will extract reads for each sample 
 * based on the `extraction.bed`, align them to the masked reference, call
 * variants using GATK HaplotypeCaller, and then combine all samples in this
 * batch (defined by $bam_set) into a single .gvcf file.
 */
process RESCUE_CAMO_VARS_PROC {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/RESCUE_CAMO_VARS", mode: 'copy')
	
	label 'RESCUE_CAMO_VARS_PROC'

	input:
		val(extraction_bed)
        val(masked_ref_fasta)
		val(gatk_bed)
        val(bam_set)
        val(max_repeats_to_rescue)

 	output:
 		path 'camo_gvcfs/*', emit: camo_gvcf_dirs

	script:

    /*
     * params.clean_tmp_files is a global variable defining whether to remove
     * tmp files for processes that could overwhelm a file system (depending on
     * how large the run is).
     */
	"""
		bash rescue_camo_variants.sh $masked_ref_fasta $extraction_bed $gatk_bed $bam_set $max_repeats_to_rescue ${params.clean_tmp_files}
	"""
}

/*
 * Prepare files and inputs for the COMBINE_AND_GENOTYPE_SAMPLES
 * process.
 */
/*
process PREP_COMBINE_AND_GENOTYPE_SAMPLES_PROC {

    label 'local'

    input:
        path gvcf_dir

    output:
        path 'split_list*', emit: filtered_lists
        tuple file(), val(), file(), file(), emit: set_of_values

    script:
    """
        for dir in ${gvcf_dirs}/*
        do
            ploidy=\$(basename \$dir)
            repeat=\$((\${ploidy##*_} / 2))

            gvcf_list="\${ploidy}.gvcfs.list"
            find \$dir -name "*.g.vcf" > \$gvcf_list

            REF=$masked_ref_fasta
            awk "\$5 == \$repeat {print \$1\":\"\$2\"-\"\$3}" \$GATK_bed | \
                while read region
                do
                    echo \$region
                    echo \$REF
                    echo \$gvcf_list
                    echo \$result_dir
                done

        done
    """
}
*/


/*
 * The process to take all .gvcfs created in RESCUE_CAMO_VARS_PROC
 * and combine them into a single .gvcf (i.e., combine all samples)
 */
/*
process COMBINE_AND_GENOTYPE_SAMPLES_PROC {

    label 'COMBINE_AND_GENOTYPE_SAMPLES'

    input:
        path(realign_bed)
        path(gatk_bed)
        path(gatk_jar)
        val(ref_prefix)
        path(filtered_list)

    output:
        path '*.dark.low_depth.bed', emit: low_depth_out
        path '*.dark.low_mapq.bed', emit: low_mapq_out

    script:
    """
        bash combine_and_genotype.sh $ $ $ $ $ $
    """
}
*/

