
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
     * Get headers for all input (b|cr)am files. 
     *
     * NOTE: 'lsamtools' is samtools compiled against 'libdeflate', which we
     * have found to be much faster than standard samtools that is compiled
     * against 'zlib'.
     */
    lsamtools_view_header_proc( sample_input_files )

    /*
     * Steps will do as follows to prepare for re-alignment:
     *   1. lsamtools_collate_and_fastq_proc: collate the input (b|cr)am files
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
    /*
    *sample_input_files
    *    | lsamtools_collate_and_fastq_proc
    *    | split_fastq_proc
    *    | map { sample_name, fastq_files ->
    *        tuple( groupKey(sample_name, fastq_files.size()), fastq_files )
    *    }
    *    | transpose()
    *    | combine( lsamtools_view_header_proc.out, by: 0 )
    *    | set { realignment_inputs_ch }
    */

    if("${params.isLongRead}" == false){
        sample_input_files
            | lsamtools_collate_and_fastq_proc
            | split_fastq_proc
            | map { sample_name, fastq_files ->
                tuple( groupKey(sample_name, fastq_files.size()), fastq_files )
            }
            | transpose()
            | combine( lsamtools_view_header_proc.out, by: 0 )
            | set { realignment_inputs_ch }

	     /*
	      * Steps will do as follows for re-alignment:
	      *   1. bwa_mem_proc: Re-align the mini .fastq files using BWA MEM
	      *   2. lsamtools_csort_proc: sort the individual mini .(b|cr)am files
	      *      generated from BWA MEM
	      *   3. groupTuple: group mini .(b|cr)am files by sample name
	      *   4. map: creates tuples of sample_name and the mini .(b|cr)am files
	      *   5. lsamtools_merge_proc: merge all of the mini .(b|cr)am files
	      */
	     bwa_mem_proc( realignment_inputs_ch )
		 | lsamtools_csort_proc
		 | groupTuple()
		 | map { sample_name, bam_files ->
		     tuple( sample_name.toString(), bam_files )
		 }
		 | lsamtools_merge_proc
		 // | view()



    } else {
         sample_input_files
            | lsamtools_fastq_proc
            | split_fastq_proc
            | map { sample_name, fastq_files ->
                tuple( groupKey(sample_name, fastq_files.size()), fastq_files )
            }
            | transpose()
            | combine( lsamtools_view_header_proc.out, by: 0 )
            | set { realignment_inputs_ch }

	
	minimap2_proc( realignment_inputs_ch )
		| lsamtools_csort_proc
		| groupTuple()
		| map{ sample_name, bam_files ->
		     tuple( sample_name.toString(), bam_files )
		 }
		 | lsamtools_merge_proc
		 // | view()
    }


    emit:
        lsamtools_merge_proc.out
}

process lsamtools_view_header_proc {

    tag { "${sample_name}" }

    label 'local'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.header.txt")

    """
    lsamtools view \\
        -H \\
        -o "${sample_name}.header.txt" \\
        "${bam}"
    """
}

process lsamtools_name_sort_proc {

    tag { "${sample_name}" }

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
    lsamtools sort \\
        -@ "${additional_threads}" \\
        -m "${mem_per_thread}" \\
        -n \\
        -o "${bam.baseName}.nsorted.bam" \\
        -T "${bam.baseName}.nsorted" \\
        "${bam}"
    """
}

process lsamtools_fastq_proc {

    tag { "${sample_name}" }

    label 'lsamtools_collate_and_fastq_proc'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.fastq.gz")

    script:

    def additional_threads = task.cpus - 1
    
    """
    lsamtools fastq \\
        -@ "${additional_threads}" \\
        -O \\
        -T RG,BC \\
        --reference "${params.original_ref}" \\
        $bam | pigz -c >"${sample_name}.fastq.gz"
    """
}


process lsamtools_collate_and_fastq_proc {

    tag { "${sample_name}" }

    label 'lsamtools_collate_and_fastq_proc'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.interleaved_R1_R2.fastq.gz")

    script:

    def additional_threads = task.cpus - 1
    
    """
    lsamtools collate \\
        -@ "${additional_threads}" \\
        -O \\
        -f \\
        "${bam}" \\
        | \\
    lsamtools fastq \\
        -@ "${additional_threads}" \\
        -O \\
        -T RG,BC \\
        -0 /dev/null \\
        --reference "${params.original_ref}" \\
        -c 1 \\
        -o "${sample_name}.interleaved_R1_R2.fastq.gz" \\
        -
    """
}

process split_fastq_proc {

    tag { "${sample_name}" }

    label 'split_fastq_proc'

    input:
    tuple val(sample_name), path(fastq)

    output:
    tuple val(sample_name), path("${sample_name}.split_${/[0-9]/*5}.fastq.gz")

    """
    lines_per_read=4
    n_lines=\$(("${params.reads_per_run}" * lines_per_read))

    pigz -dcp 4 "${fastq}" \\
    | \\
    split \\
    --suffix-length=5 \\
    --additional-suffix=".fastq" \\
    --filter='pigz --fast -p 12 > \$FILE.gz' \\
    -d \\
    -l \$n_lines \\
    - \\
    "${sample_name}.split_"
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
    if [ "${params.output_format}" = "bam" ]; then
        bwa mem \\
            -p \\
            -t "${task_cpus}" \\
            -M \\
            -C \\
            -H <(grep "^@RG" "${header}") \\
            "${params.align_to_ref}" \\
            "${fastq}" |
        lsamtools view \\
            -1 \\
            -o "${fastq.baseName}.unsorted.mini.bam" \\
            -
    else
        bwa mem \\
            -p \\
            -t "${task_cpus}" \\
            -M \\
            -C \\
            -H <(grep "^@RG" "${header}") \\
            "${params.align_to_ref}" \\
            "${fastq}" |
        lsamtools view \\
            -C \\
            -T params.align_to_ref \\
            -o "${fastq.baseName}.unsorted.mini.cram" \\
            -
    fi
    """
}


process minimap2_proc {

    tag { "${sample_name}:${fastq.name}" }

    label 'minimap2_proc'

    input:
    tuple val(sample_name), path(fastq), path(header)

    output:
    tuple val(sample_name), path("${fastq.baseName}.unsorted.mini.bam")

    script:
    def task_cpus = task.cpus > 1 ? task.cpus - 1 : task.cpus

    print sample_name
    print fastq
    """
   
    if [ "${params.isONT}" = "true" ]; then
        if [ "${params.output_format}" = "bam" ]; then
            minimap2 \\
                -t "${task_cpus}" \\
                -a \\
                -x map-ont \\
                -Y \\
                --eqx \\
                --secondary=no \\
                -L \\
                -O 5,56 \\
                -E 4,1 \\
                -B 5 \\
                -z 400,50 \\
                -r 2k \\
                "${params.align_to_ref}" \\
                "${fastq}" |
            lsamtools view \\
                -1 \\
                -O bam \\
                -o "${fastq.baseName}.unsorted.mini.bam" \\
                -
        else
            minimap2 \\
                -t "${task_cpus}" \\
                -a \\
                -Y \\
                --eqx \\
                --secondary=no \\
                -L \\
                -O 5,56 \\
                -E 4,1 \\
                -B 5 \\
                -z 400,50 \\
                -r 2k \\
                -x map-ont \\
                "${params.align_to_ref}" \\
                "${fastq}" |
            lsamtools view \\
                -C \\
                -T "${params.align_to_ref}" \\
                -o "${fastq.baseName}.unsorted.mini.cram" \\
                -
        fi
    elif [ "${params.isPacBio}" = "true" ]; then
        if [ "${params.output_format}" = "bam" ]; then
            minimap2 \\
                -t "${task_cpus}" \\
                -a \\
                -x map-pb \\
                -Y \\
                --eqx \\
                --secondary=no \\
                -L \\
                -O 5,56 \\
                -E 4,1 \\
                -B 5 \\
                -z 400,50 \\
                -r 2k \\
                "${params.align_to_ref}" \\
                "${fastq}" |
            lsamtools view \\
                -1 \\
                -O bam \\
                -o "${fastq.baseName}.unsorted.mini.bam" \\
                -
        else
            minimap2 \\
                -t "${task_cpus}" \\
                -a \\
                -x map-pb \\ 
                --secondary=no \\
                -Y \\
                --eqx \\
                -L \\
                -O 5,56 \\
                -E 4,1 \\
                -B 5 \\
                -z 400,50 \\
                -r 2k \\
                "${params.align_to_ref}" \\
                "${fastq}" |
            lsamtools view \\
                -C \\
                -T "${params.align_to_ref}" \\
                -o "${fastq.baseName}.unsorted.mini.cram" \\
                -
        fi
    fi
    """
}

process csort_proc {

    tag { "${sample_name}:${bam.name}" }

    label 'csort_proc'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${bam.baseName}.csorted.{bam,cram}")

    script:

    def additional_threads = task.cpus - 1

    /*
     * Calculate mem per thread. Multiply total mem by 0.8 to provide a 20%
     * buffer.
     */
    def avail_mem = task.memory ? ( task.memory.toGiga() * 0.8 ).toInteger().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "${avail_mem}G" : ''

    """
    echo "Available mem: ${avail_mem}"
    echo "Mem per thread: ${mem_per_thread}"

    if [ "${params.output_format}" = "bam" ]; then
        sambamba sort \\
            -t "${task.cpus}" \\
            -f bam \\
            -m "${avail_mem}" \\
            -o "${bam.baseName}.csorted.bam" \\
            -l 2 \\
            "${bam}"

        sambamba index \\
            -t "${task.cpus}" \\
            --show-progress \\
            "${bam.baseName}.csorted.bam"

    else
        # sambamba does not support .cram, so still using lsamtools
        lsamtools sort \\
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

process lsamtools_csort_proc {

    tag { "${sample_name}:${bam.name}" }

    label 'lsamtools_csort_proc'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${bam.baseName}.csorted.{bam,cram}")

    script:

    def additional_threads = task.cpus - 1

    /*
     * Calculate mem per thread. Multiply total mem by 0.8 to provide a 20%
     * buffer. Multiply by 1024 to convert to MB
     */
    def avail_mem = task.memory ? ( task.memory.toGiga() * 0.8 * 1024 ).toInteger().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "${avail_mem}M" : ''

    """
    echo "Available mem: $avail_mem"
    echo "Mem per thread: $mem_per_thread"
    if [ "${params.output_format}" = "bam" ]; then
        lsamtools sort \\
            -@ "${additional_threads}" \\
            -m "${mem_per_thread}" \\
            -o "${bam.baseName}.csorted.bam" \\
            -T "${bam.baseName}.csorted" \\
            --write-index \\
            "${bam}"
    else
        lsamtools sort \\
            -@ "${additional_threads}" \\
            -m "${mem_per_thread}" \\
            -o "${bam.baseName}.csorted.bam" \\
            -T "${bam.baseName}.csorted" \\
            -O cram \\
            --reference "${params.align_to_ref}" \\
            --write-index \\
            "${bam}"
    fi
    """
}


process lsamtools_merge_proc {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/01-REALIGN_SAMPLES", mode: 'copy')

    tag { sample_name }

    label 'lsamtools_merge_proc'

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name),
            path("${sample_name}*.sorted.merged.final.{bam,cram}"),
            path("${sample_name}*.sorted.merged.final.{bam,cram}{.bai,.csi,.crai}")

    script:

    def additional_threads = task.cpus - 1

    """
    if [ "${params.output_format}" = "bam" ]; then
        lsamtools merge \\
            -@ "${additional_threads}" \\
            -c \\
            -p \\
            "${sample_name}.${params.align_to_ref_tag}.sorted.merged.final.bam" \\
            --write-index \\
            ${bams}
    else
        lsamtools merge \\
            -@ "${additional_threads}" \\
            -c \\
            -p \\
            -O cram \\
            "${sample_name}.${params.align_to_ref_tag}.sorted.merged.final.cram" \\
            --write-index \\
            ${bams}
    fi
    """
}

