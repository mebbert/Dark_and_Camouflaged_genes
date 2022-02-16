
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2


/*
 * The code herein is a complete conversion of our original .bash scripts to
 * Nextflow. It was converted by @Steve on bioinformatics.stackexchange.com:
 * https://bioinformatics.stackexchange.com/questions/18549/nextflow-dsl-v2-how-to-best-synchronize-multiple-outputs-from-a-single-proces/18566#18566
 *
 * I (Mark Ebbert) have modified the code a bit to handle .cram files and change
 * naming conventions.
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

    samtools_view_header_proc( sample_input_files )

    sample_input_files
        | samtools_collate_and_fastq_proc
        | split_fastq_proc
        | map { sample_name, fastq_files ->
            tuple( groupKey(sample_name, fastq_files.size()), fastq_files )
        }
        | view()
        | transpose()
        | view()
        | combine( samtools_view_header_proc.out, by: 0 )
        | view()
        | set { realignment_inputs }

     bwa_mem_proc( realignment_inputs )
         | samtools_coordinate_sort_proc
         | groupTuple()
         | map { sample_name, bam_files ->
             tuple( sample_name.toString(), bam_files )
         }
         | samtools_merge_proc
         | view()
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
    def avail_mem = task.memory ? task.memory.toGiga().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "-m ${avail_mem}G" : ''

    """
    samtools sort \\
        -@ "${task.cpus - 1}" \\
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

    """
    samtools collate \\
        -@ "${task.cpus - 1}" \\
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
    split_fastq.py $sample_name $fastq $params.reads_per_run
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
    tuple val(sample_name), path("${bam.baseName}.csorted.bam")

    script:
    def avail_mem = task.memory ? task.memory.toGiga().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "${avail_mem}G" : ''

    """
    if [ "$params.output_format" = "bam" ]; then
        samtools sort \\
            -@ "${task.cpus - 1}" \\
            -m ${mem_per_thread} \\
            -o "${bam.baseName}.csorted.bam" \\
            -T "${bam.baseName}.csorted" \\
            --write-index \\
            "${bam}"
    else
        samtools sort \\
            -@ "${task.cpus - 1}" \\
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

    tag { sample_name }

    label 'samtools_merge_proc'

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name), path("${sample_name}.bam{,.bai}")

    """
    samtools merge \\
        -c \\
        -p \\
        "${sample_name}.bam" \\
        ${bams}
    samtools index \\
        "${sample_name}.bam"
    """
}

