import groovy.io.FileType

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

params.GVCF_DIRECTORY = './path/to/directories'

workflow RESCUE_CAMO_VARS_WF {
    take:
        bam_path
        extraction_bed
        masked_ref_fasta
        gatk_bed
        n_samples_per_batch
        max_repeats_to_rescue

    main:

        PREP_BATCH_FILES_PROC(bam_path, n_samples_per_batch)
        RESCUE_CAMO_VARS_PROC(
                                extraction_bed,
                                masked_ref_fasta,
                                gatk_bed,
                                PREP_BATCH_FILES_PROC.out.bam_sets.flatten(),
                                max_repeats_to_rescue
                              )

    /*
     * Prepare input parameters for COMBINE_AND_GENOTYPE
     */
    Channel.fromPath( gatk_bed ) \
        | splitCsv(sep: '\t') \
        | map { line ->
            tuple( line[4].toInteger(), "${line[0]}:${line[1]}-${line[2]}" )
        } \
        | set { regions }

//     test = ""
// 
//     RESCUE_CAMO_VARS_PROC.out.camo_gvcf_parent_dir.concatenate() // /**/*.g.vcf"
//     RESCUE_CAMO_VARS_PROC.out.camo_gvcf_parent_dir.subscribe{ test += it.name } // /**/*.g.vcf"
// 
//     dfile = file('debug.txt')
//     dfile.append("GVCF parent dir: \n")
//     dfile.append(RESCUE_CAMO_VARS_PROC.out.camo_gvcf_parent_dir.view().toString() + "\n")
//     dfile.append(params.GVCF_DIRECTORY.toString() + "\n")
//     dfile.append("Test: " + test + "\n")
//     dfile.append(test)
    
    /*
     * TODO: Figure out how to access the emitted path to the .gvcfs
     * from RESCUE_CAMO_VARS_PROC so that we don't have to rely on the
     * saved results and hard code the directory.
     */
    GVCF_dir = file("${params.results_dir}/RESCUE_CAMO_VARS")
    c_and_g_input_params = Channel.fromPath( "${GVCF_dir}/**/*.g.vcf" ) \
         .map { tuple( GVCF_dir.relativize(it).subpath(1,2).name, it ) }
         .groupTuple()
         .map { dirname, gvcf_files ->
             def ploidy = dirname.replaceFirst(/^.*_/, "").toInteger()
             def repeat = ploidy.intdiv(2)
             def ref_fasta = file( masked_ref_fasta )
             tuple( repeat, dirname, ref_fasta, gvcf_files )
         }
         .combine( regions, by: 0 )
         .map { repeat, dirname, ref_fasta, gvcf_files, region ->
             tuple( dirname, region, ref_fasta, gvcf_files )
         }
         
    // c_and_g_input_params.view()
    // COMBINE_AND_GENOTYPE_PROC(c_and_g_input_params, RESCUE_CAMO_VARS_PROC.out.rescue_complete)
    COMBINE_AND_GENOTYPE_PROC(c_and_g_input_params)
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
 		path 'camo_batch_gvcfs/', emit: camo_gvcf_parent_dir
        val 'complete', emit: rescue_complete

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
 *
 */
process COMBINE_AND_GENOTYPE_PROC {

    publishDir "${params.results_dir}/COMBINE_AND_GENOTYPE/${ploidy_group}"

	label 'COMBINE_AND_GENOTYPE'

    input:

    /*
     * Receiving as ref_fasta and gvcf_files as 'val' input so that NextFlow
     * will not 'stage' the files in the working directory because it doesn't
     * stage the index files with them. The other option would be to also receive
     * the index files as paths so both are staged.
     */
    tuple val(ploidy_group), val(region_string), val(ref_fasta), val(gvcf_files)

    /*
     * Passing in this value from RESCUE_CAMO_VARS_PROC to signal that
     * this process must wait for the rescue to complete.
     */
    // val(rescue_complete) 

    output:
    tuple val(ploidy_group), val(region_string), path("full_cohort.combined.${region}.g.vcf")

    script:

    region = region_string.replaceAll(':', '_')

    def avail_mem = task.memory ? task.memory.toGiga() : 0

    def Xmx = avail_mem >= 8 ? "-Xmx${avail_mem - 1}G" : ''
    def Xms = avail_mem >= 8 ? "-Xms${avail_mem.intdiv(2)}G" : ''

    """

    echo "Ploidy group: $ploidy_group"
    echo "Region: $region_string"
    echo "Ref: $ref_fasta"

    for gvcf_file in $gvcf_files
    do
        echo "GVCF file: \$gvcf_file"
    done

    cat << __EOF__ > "${ploidy_group}.gvcf.list"
    ${gvcf_files.join('\n'+' '*4)}
    __EOF__

    gatk \\
        --java-options "${Xmx} ${Xms} -XX:+UseSerialGC" \\
        CombineGVCFs \\
        -R "${ref_fasta}" \\
        -L "${region_string}" \\
        -O "full_cohort.combined.${region}.g.vcf" \\
        -V "${ploidy_group}.gvcf.list"

    gatk \\
        --java-options "${Xmx} ${Xms} -XX:+UseSerialGC" \\
        GenotypeGVCFs \\
        -R "${ref_fasta}" \\
        -L "${region_string}" \\
        -O "full_cohort.combined.${region}.vcf" \\
        -V "full_cohort.combined.${region}.g.vcf" \\
        -A GenotypeSummaries
    """
}

