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
        ref_fasta
        masked_ref_fasta
        gatk_bed
        n_samples_per_batch
        max_repeats_to_rescue
        camo_annotations
        align_to_bed
        ref_tag

    main:

        /*
         * Prep batch files for RESCUE_CAMO_VARS_PROC
         */
        PREP_BATCH_FILES_PROC(bam_path, n_samples_per_batch)

        /*
         * Rescue camo vars
         */
        RESCUE_CAMO_VARS_PROC(
                                extraction_bed,
                                masked_ref_fasta,
                                gatk_bed,
                                PREP_BATCH_FILES_PROC.out.bam_sets.flatten(),
                                max_repeats_to_rescue
                              )

        /*
         * Prepare input parameters for COMBINE_AND_GENOTYPE and run. Most of this 
         * `NextFlow/Groovy` wizardry is courtesy of @Steve on `StackOverflow`
         * https://stackoverflow.com/questions/70718115/the-most-nextflow-like-dsl2-way-to-incorporate-a-former-bash-scheduler-submiss/70735565#70735565
         */
        Channel.fromPath( gatk_bed ) \
            | splitCsv(sep: '\t') \
            | map { line ->
                tuple( line[4].toInteger(), "${line[0]}:${line[1]}-${line[2]}" )
            } \
            | set { regions }

        RESCUE_CAMO_VARS_PROC.out.camo_gvcfs.flatten()
             .map { tuple( it.subpath(it.getNameCount() - 2, it.getNameCount() - 1).name, it )
             }
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
             } | COMBINE_AND_GENOTYPE_PROC


        /*
         * Identify false positives (i.e., reference-based artifiacts).
         */
        IDENTIFY_FALSE_POSITIVES_PROC(camo_annotations, align_to_bed, ref_fasta, ref_tag)
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
 * The process to rescue camouflaged variants. This process performs the
 * following steps (for each batch of samples):
 *   1. Extract low MAPQ reads (MAPQ < 10) for each sample based on the
 *      `extraction.bed`
 *   2. Align extracted reads to the masked reference created in MASK_GENOME
 *   3. Call variants using GATK HaplotypeCaller with the GATK .bed file
 *   4. Combine all samples in this batch (defined by $bam_set) and for each
 *      ploidy into a single .gvcf file.
 *
 * This is done seperately for each repeat number, adjusting `ploidy` for GATK
 * HaplotypeCaller for the specific region being called (i.e. regions with
 * repeat number 2 uses ploidy of 4). All samples get combined into ploidy-
 * specific GVCFs for each ploidy tested.
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
        /*
         * Emit the full path to all .gvcf files
         */
 		// path 'camo_batch_gvcfs/', emit: camo_gvcfs
 		path 'camo_batch_gvcfs/**/*.g.vcf', emit: camo_gvcfs
        // val 'complete', emit: rescue_complete

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
 * The process to combine and genotype all samples for each ploidy .gvcf
 * generated in `RESCUE_CAMO_VARS_PROC`. Combining all gVCFs across all camo
 * regions for a large cohort (e.g., the ADSP) requires too much memory, so
 * the job is parallelized by looking at each camo region, individually, and
 * combine all samples and genotype them for just this single region.
 *
 * In a later step, variants from all camo region VCFs are catted into a
 * combined cohort VCF for a given ploidy across all samples.
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

    combined_region_gvcf="full_cohort.combined.${region}.g.vcf"
    final_region_vcf="full_cohort.combined.${region}.vcf"
    input_gvcf_list="${ploidy_group}.gvcf.list"

    gatk \\
        --java-options "${Xmx} ${Xms} -XX:+UseSerialGC" \\
        CombineGVCFs \\
        -R "${ref_fasta}" \\
        -L "${region_string}" \\
        -O \$combined_region_gvcf \\
        -V \$input_gvcf_list

    gatk \\
        --java-options "${Xmx} ${Xms} -XX:+UseSerialGC" \\
        GenotypeGVCFs \\
        -R "${ref_fasta}" \\
        -L "${region_string}" \\
        -A GenotypeSummaries \\
        -O \$final_region_vcf \\
        -V $combined_region_gvcf
    """
}


/*
 * The process to identify false-positive variant calls in the rescued camo
 * variants. These false-positives specifically originate from what we term
 * "reference-based artifacts". These artifacts show up when when two regions
 * in a given camo set are not 100% identical. Thus, when a read from one
 * region is forced to align to the other, an artificial variant is called.
 *
 * This process creates a list of all possible reference-based artifacts by
 * taking all camo CDS regions with repeat number ≤ 5 (i.e., ≤ 5 camo regions
 * in the set) and blat them against the whole genome. Any DNA sequence from a
 * hit with sequence identity ≥ 98% are locally realigned back to the query
 * seqeunce using Bio.pairwise2. Any mismatches or gaps in the aligned sequence
 * are converted into variant positions and output to a bed file that lists
 * positions of reference-based-artifacts to be filtered from the VCF files.
 *
 * This only needs to be run once per reference genome.
 */
process IDENTIFY_FALSE_POSITIVES_PROC {

    publishDir "${params.results_dir}/IDENTIFY_FALSE_POSITIVES"

	label 'IDENIFY_FALSE_POSITIVES'

    input:
        path camo_annotations
        path align_to_bed
        val ref_fasta
        val ref_name
        
    output:
        path "reference_based_artifacts.${ref_name}.bed"

    script:
    """
    artifacts_bed="reference_based_artifacts.${ref_name}.bed"

    bash identify_false_positives.sh $camo_annotations $align_to_bed $ref_fasta \$artifacts_bed $task.cpus
    """
}
