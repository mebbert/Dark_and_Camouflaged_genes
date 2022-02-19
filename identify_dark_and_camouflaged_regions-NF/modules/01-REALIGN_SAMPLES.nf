
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2


/*
 * The code herein is a complete conversion of our original .bash scripts to
 * Nextflow. It was converted by @Steve on bioinformatics.stackexchange.com:
 * https://bioinformatics.stackexchange.com/questions/18549/nextflow-dsl-v2-how-to-best-synchronize-multiple-outputs-from-a-single-proces/18566#18566
 *
 * @Steve has been a huge help. I (Mark Ebbert) modified the code a bit to handle
 * .cram files and change naming conventions.
 */


/*
 * REALIGN_SAMPLES_WF is the main workflow for re-aligning .(cr|b)am files
 */
workflow REALIGN_SAMPLES_WF {

    take:
        sample_input_files_path


    main:

    /*
     * Create sample input tuples
     */
    Channel.fromPath(sample_input_files_path, checkIfExists: true)
        | filter( ~/.*(\.sam|\.bam|\.cram)/ )
        | map { tuple( it.baseName, it ) }
        | set { sample_input_files }

    /*
     * Get headers for all input (b|cr)am files
     */
    samtools_view_header_proc( sample_input_files )

    /*
     * Steps will do as follows to prepare for re-alignment:
     *   1. samtools_collate_and_fastq_proc: collate the input (b|cr)am files
     *      and convert to a single, interleaved .fastq file
     *   2. split_fastq_proc: split the single .fastq file into chunks of
     *      'params.reads_per_run' reads for BWA alignment
     *   3. map: create tuples of sample_name, the number of .fastqs generated
     *      for the given sample, and the actual .fastq files.
     *   4. transpose: transposes the set of tuples
     *   5. combine: adds the header files to the tuples
     *   6. set: creates an output channel of the tuples
     *
     *
     * NOTE: If 'params.reads_per_run' results in a single .fastq file for a
     * given sample, Nextflow will error out with something like the following:
     * "Invalid method invocation `groupKey` with arguments: <sample_name>
     *    (java.lang.String), 26766997 (java.lang.Long) on Nextflow type"
     *
     * This error happens because `fastq_files.size()` returns the file size for the single file
     * returned rather than reporting the size of the list of files. We expect.
     * a *list* of files. Nextflow channels seem inherently buggy to me.
     *
     * To 'fix' it, reduce params.reads_per_run to a number less than the number
     * of paired-end reads
     */
    sample_input_files
        | samtools_collate_and_fastq_proc
        | split_fastq_proc
        | map { sample_name, fastq_files ->
            tuple( groupKey(sample_name, fastq_files.size()), fastq_files )
        }
        | transpose()
        | combine( samtools_view_header_proc.out, by: 0 )
        | set { realignment_inputs_ch }


     /*
      * Steps will do as follows for re-alignment:
      *   1. bwa_mem_proc: Re-align the mini .fastq files using BWA MEM
      *   2. samtools_coordinate_sort_proc: sort the individual mini .(b|cr)am files
      *      generated from BWA MEM
      *   3. groupTuple: group mini .(b|cr)am files by sample name
      *   4. map: creates tuples of sample_name and the mini .(b|cr)am files
      *   5. samtools_merge_proc: merge all of the mini .(b|cr)am files
      */
     bwa_mem_proc( realignment_inputs_ch )
         | samtools_coordinate_sort_proc
         | groupTuple()
         | map { sample_name, bam_files ->
             tuple( sample_name.toString(), bam_files )
         }
         | samtools_merge_proc
         | view()

    emit:
        samtools_merge_proc.out
}

process samtools_view_header_proc {

    tag { "${sample_name}:${bam.name}" }

    label 'local'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.header.txt")

    """
    samtools view \\
        -H \\
        -o "${sample_name}.header.txt" \\
        "${bam}"
    """
}

process samtools_name_sort_proc {

    tag { "${sample_name}:${bam.name}" }

    cpus 4
    memory 8.GB

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${bam.baseName}.nsorted.bam")

    script:

    def additional_threads = task.cpus - 1

    def avail_mem = task.memory ? task.memory.toGiga().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "-m ${avail_mem}G" : ''

    """
    samtools sort \\
        -@ "${additional_threads}" \\
        ${mem_per_thread} \\
        -n \\
        -o "${bam.baseName}.nsorted.bam" \\
        -T "${bam.baseName}.nsorted" \\
        "${bam}"
    """
}

process samtools_collate_and_fastq_proc {

    tag { "${sample_name}:${bam.name}" }

    label 'samtools_collate_and_fastq_proc'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.interleaved_R1_R2.fastq.gz")

    script:

    def additional_threads = task.cpus - 1
    
    """
    samtools collate \\
        -@ "${additional_threads}" \\
        -O \\
        -f \\
        "${bam}" \\
        | \\
    samtools fastq \\
        -O \\
        -T RG,BC \\
        -0 /dev/null \\
        --reference $params.original_ref \\
        -c 1 \\
        -o "${sample_name}.interleaved_R1_R2.fastq.gz" \\
        -
    """
}

process split_fastq_proc {

    tag { "${sample_name}:${fastq.name}" }

    label 'split_fastq_proc'

    input:
    tuple val(sample_name), path(fastq)

    output:
    tuple val(sample_name), path("${sample_name}.interleaved_R1_R2.split_${/[0-9]/*5}.fastq.gz")

    """
    lines_per_read=4

    zcat "${fastq}" \\
    | \\
    split \\
    --suffix-length=5 \\
    --additional-suffix=".fastq" \\
    --filter='gzip > \$FILE.gz' \\
    -d \\
    -l "${params.reads_per_run * lines_per_read}" \\
    - \\
    "${sample_name}.interleaved_R1_R2.split_"
    """
}

process bwa_index_proc {

    tag { fasta.name }

    cpus 1
    memory 12.GB

    input:
    path fasta

    output:
    tuple val(fasta.name), path("${fasta}.{amb,ann,bwt,pac,sa}")

    """
    bwa index "${fasta}"
    """
}

process bwa_mem_proc {

    tag { "${sample_name}:${fastq.name}" }

    label 'bwa_mem_proc'

    input:
    tuple val(sample_name), path(fastq), path(header)

    output:
    tuple val(sample_name), path("${fastq.baseName}.unsorted.mini.bam")

    script:
    def task_cpus = task.cpus > 1 ? task.cpus - 1 : task.cpus

    """
    if [ "$params.output_format" = "bam" ]; then
        bwa mem \\
            -p \\
            -t ${task_cpus} \\
            -M \\
            -C \\
            -H <(grep "^@RG" "${header}") \\
            "${params.align_to_ref}" \\
            "${fastq}" |
        samtools view \\
            -1 \\
            -o "${fastq.baseName}.unsorted.mini.bam" \\
            -
    else
        bwa mem \\
            -p \\
            -t ${task_cpus} \\
            -M \\
            -C \\
            -H <(grep "^@RG" "${header}") \\
            "${params.align_to_ref}" \\
            "${fastq}" |
        samtools view \\
            -C \\
            -T params.align_to_ref \\
            -o "${fastq.baseName}.unsorted.mini.bam" \\
            -
    fi
    """
}

process samtools_coordinate_sort_proc {

    tag { "${sample_name}:${bam.name}" }

    label 'samtools_coordinate_sort_proc'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${bam.baseName}.csorted.{bam,cram}")

    script:

    def additional_threads = task.cpus - 1

    def avail_mem = task.memory ? task.memory.toGiga().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "${avail_mem}G" : ''

    """
    if [ "$params.output_format" = "bam" ]; then
        samtools sort \\
            -@ "${additional_threads}" \\
            -m ${mem_per_thread} \\
            -o "${bam.baseName}.csorted.bam" \\
            -T "${bam.baseName}.csorted" \\
            --write-index \\
            "${bam}"
    else
        samtools sort \\
            -@ "${additional_threads}" \\
            -m ${mem_per_thread} \\
            -o "${bam.baseName}.csorted.bam" \\
            -T "${bam.baseName}.csorted" \\
            -O cram \\
            --reference "${params.align_to_ref}" \\
            --write-index \\
            "${bam}"
    fi
    """
}

process samtools_merge_proc {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/01-REALIGN_SAMPLES", mode: 'copy')

    tag { sample_name }

    label 'samtools_merge_proc'

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name), path("${sample_name}*.sorted.merged.final.{bam,cram}{,.bai}")

    """
    if [ "$params.output_format" = "bam" ]; then
        samtools merge \\
            -c \\
            -p \\
            "${sample_name}.${params.align_to_ref_tag}.sorted.merged.final.bam" \\
            --write-index \\
            ${bams}
    else
        samtools merge \\
            -c \\
            -p \\
            -O cram \\
            "${sample_name}.${params.align_to_ref_tag}.sorted.merged.final.cram" \\
            --write-index \\
            ${bams}
    fi
    """
}

