
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

workflow REALIGN_SAMPLES_WF {

    take:
        sample_input_file_ch
        n_reads_per_run
        original_ref
        align_to_ref
        align_to_ref_tag
        output_format

    main:

    /*
     * Generate fastq files from the original .(cr|b)am file and split into
     * mini .fastq files with `n_reads_per_run` in each mini .fastq
     */

    // sample_input_file_ch.view()
    GENERATE_FASTQS_PROC(sample_input_file_ch, n_reads_per_run, original_ref)


    /*
     * Get sample name and read group from the RG file using Groovy
     */

     GENERATE_FASTQS_PROC.out.sample_tuple_file
        .map{
            sample_tuple_file ->
                sample_tuple_file.withReader {
                    it.eachLine {
                        line ->
                            println "line: " + line
                            toks = line.split(',')
                            sample_name = toks[0]
                            fq = toks[1]
                            rg = toks[2]

                            tuple(sample_name, fq, rg)
                    }
                }
        }
        .view()
        .set { sample_fq_tuples_ch }



//    sample_name_rg = GENERATE_FASTQS_PROC.out.sample_RG_file
//                    .map {
//                        rg_file ->
//                            file_lines = rg_file.readLines()
//                            sample_name = file_lines[0]
//                            sample_RG = file_lines[1]
//
//                            tuple( sample_name, sample_RG )
//                    }


//    sample_name_and_RG_ch = GENERATE_FASTQS_PROC.out.rg_file_path
//                                .map{
//
//                                    rg_file ->
//                                        
//                                        /* get sample name and @RG from file */
//                                        file_lines = rg_file.readLines()
//
//                                        sample_name = file_lines[0]
//                                        sample_RG = file_lines[1]
//
//                                        [ sample_name, sample_RG ]
//                                }
//                                .view()


//    fq_split_reads = GENERATE_FASTQS_PROC.out.fastq_path
//                        .view()
//                        .map {
//                            split_fq_path ->
//                                println "split_fq_path: " + split_fq_path
//                                Channel.fromPath(split_fq_path.toString())
//                                .subscribe { print it.flatten() }
//                                        fq_num = fq.getName().split('\\.')[-2]
//                                        tuple( fq_num, fq )
//                        }
//                        .flatten()
//                        .view()

//     GENERATE_FASTQS_PROC.out.fastqs
//        .flatten()
//        .combine(sample_name_and_RG_ch)
//        .view()
//        .set{fq_split_reads_and_sample_name_RG}

    /*
     * Group the split .fastq files into appropriate pairs and combine them
     * with the sample Read Group (RG) from the original .(cr|b)am file.
     *
     * There's probably a simpler way to do this, but the level of abstraction
     * in Nextflow (combined with insufficient documentation) is driving me nuts.
     *
     * Here are the steps explained:
     *  1.   .flatten: get individual elements
     *  2.       .map: loop over each element left in the channel (i.e., the .fastq
     *                 files), extract the number assigned to the file pair by
     *                 `split`, and make a new tuple. The 'fq_num' is used to 
     *                 pair the files.
     *  3. .groupTuple: group the file pairs by `fq_num`
     *  4.   .combine: add the sample RG to the tuple based on the
     *                 `sample_input_file.baseName` ID
     *  5.       .map: the .combine leaves a tuple of tuples. For each individual
     *                 tuple, flatten and create into a 'clean' tuple
     */
//    fq_split_reads_and_RG_ch = GENERATE_FASTQS_PROC.out.RG_and_fastqs
//                                .flatten()
//                                .view()
//                                .map {
//                                    file_path ->
//                                        println item
//                                        file_name = file_path.getName()
//                                        if(file_name.endsWith('*.txt')){
//                                            file_lines = rg_file.readLines()
//                                            sample_name = file_lines[0]
//                                            sample_RG = file_lines[1]
//                                        }
//
//                                        tuple(fq1, fq2, rg_file)
//                                        
//                                }
//                                .flatten()
//                                .map {
//                                    fq ->
//
//                                        /*
//                                         * Get the second to last item--the matching number from
//                                         * split
//                                         */
//                                        fq_num = fq.getName().split('\\.')[-2]
//                                        tuple( fq_num, fq )
//                                }
//                                .groupTuple()
//                                .view()
//                                .combine( sample_name_rg )
//                                .view()
//                                .map {
//                                    input_tuple ->
//                                        tuple( input_tuple.flatten() )
//                                }
//                                .view()


    /*
     * Align fastq file pairs
     */

    ALIGN_FASTQ_BWA_PROC(sample_fq_tuples_ch, align_to_ref,
                             align_to_ref_tag, original_ref, output_format)
//
//
//
//    /*
//     * Sort and index the mini .(cr|b)am files
//     */
//    SAMTOOLS_SORT_AND_INDEX_MINI_PROC(ALIGN_FASTQ_BWA_PROC.out.final_alignment,
//                                ALIGN_FASTQ_BWA_PROC.out.sample_name, output_format)
//
//
//    /*
//     * Merge the mini .(cr|b)am files
//     */
//    /* Collect all output from RUN_DRF_PROC and group by sample */
//    samtools_sort_output_by_sample = SAMTOOLS_SORT_AND_INDEX_MINI_PROC.out.sample_name_and_mini_output.collect().flatten()
//                    .map {
//                        sample_file_path ->
//                            // println "sample file path: " + sample_file_path
//                            // println "sample file path name: " + sample_file_path.name
//                            sample_name = sample_file_path.name.split('\\.')[0]
//                            tuple(
//                                    sample_name, sample_file_path
//                                 )
//                    }
//                    .groupTuple()
//
//
//    SAMTOOLS_MERGE_AND_INDEX_PROC(samtools_sort_output_by_sample, align_to_ref_tag, output_format)

}


process GENERATE_FASTQS_PROC {

    label 'GENERATE_FASTQS_PROC'

    input:
        path(sample_input_file)
        val(n_reads_per_run)
        val(original_ref)

    output:
        path("*.tuples.txt"), emit: sample_tuple_file
        // path("*.ReadGroup.txt"), emit: rg_file_path
        // path("split_fqs/*.fastq"), emit: fastqs
        // tuple path("*.ReadGroup.txt"), path("split_fqs/*.fastq"), \
        //    emit: RG_and_fastqs
        // tuple val(sample_input_file.baseName), path("*.ReadGroup.txt"), \
        //    emit: sample_RG_file

    script:

    /*
     * Calculate the mem per thread. Divide by double the number of allocated threads to provide
     * buffer. Samtools' memory limit is wildly inaccurate and unpredictable, in my experience.
     */
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def mem_per_thread = ( avail_mem ).intdiv( task.cpus * 2 )


    """
    generate_fastqs.sh $sample_input_file $n_reads_per_run $original_ref $task.cpus $mem_per_thread
    """

}



process ALIGN_FASTQ_BWA_PROC {

	label 'REALIGN_SAMPLES'

	input:
        // tuple val(split_ID), path(fq1), path(fq2), path(sample_RG_file)
        // tuple val(split_ID), path(fq1), path(fq2), val(sample_name), val(sampe_RG)
        tuple val(sample_name), path(fq), val(sample_RG)
		val(align_to_ref)
		val(align_to_ref_tag)
		val(original_ref)
        val(output_format)


	output:
        path('*unsorted.mini.{bam,cram}'), emit: final_alignment
        val(sample_name)

	script:

    /*
     * Sanity check: verify the split_ID (provided by Unix `split` when splitting
     * the files) matches the split ID within the file names.
     *
     * Also verify provided sample name matches file names
     */

//     fq1_name = fq1.getName()
//     fq2_name = fq2.getName()
//     fq1_split_ID = fq1_name.split('\\.')[-2].toInteger()
//     fq2_split_ID = fq2_name.split('\\.')[-2].toInteger()
// 
//     match = fq1_name =~ /([A-Za-z0-9_\-]+)_R[12]\.split\.\d+.fastq/
//     fq1_sample_name = match[0][1]
// 
//     match = fq2_name =~ /([A-Za-z0-9_\-]+)_R[12]\.split\.\d+.fastq/
//     fq2_sample_name = match[0][1]
// 
//     println fq1_sample_name
//     println fq2_sample_name
// 
//     if( split_ID.toInteger() != fq1_split_ID ||
//             split_ID.toInteger() != fq2_split_ID ){
//         
//         throw new Exception("ERROR: The .fastq files received in ALIGN_FASTQ_BWA_PROC do not" +
//                 " match the split_ID provided.\n" +
//                 "Split_ID received: $split_ID\n" +
//                 "Fastq 1: $fq1\n" +
//                 "Fastq 2: $fq2\n"
//                 )
//     }
// 
//     if( !fq1_sample_name.equals(sample_name) ||
//             !fq2_sample_name.equals(sample_name) ){
// 
//         throw new Exception("ERROR: The .fastq files received in ALIGN_FASTQ_BWA_PROC do not" +
//                 " match the sample name provided.\n" +
//                 "Sample name received: $sample_name\n" +
//                 "Fastq 1: $fq1\n" +
//                 "Fastq 2: $fq2\n"
//                 )
//     }

	"""
   #  lines=()
   #  while read line
   #  do
   #      lines+=("\$line")
   #  done < $sample_RG_file

   #  sample_name="\${lines[0]}"
   #  sample_RG="\${lines[1]}"

   #  echo "Sample RG: \$sample_RG"

   #  echo "Sample fq1: $fq1"
   #  echo "Sample fq2: $fq2"
	bash realign_bwa-split.sh $fq \$sample_name \"\$sample_RG\" ${align_to_ref} ${output_format} $task.cpus
	"""
}




process SAMTOOLS_SORT_AND_INDEX_MINI_PROC {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    // publishDir("${params.results_dir}/01-REALIGN_SAMPLES", mode: 'copy')

    label 'SORT_AND_INDEX_SAMPLE_PROC'

    input:
        path(final_mini_alignment_file)
        val(sample_name)
        val(output_format)

    output:
        tuple(val(sample_name), path '*sorted.mini.{bam,cram}*'), emit: sample_name_and_mini_output

    script:

    /*
     * Calculate the mem per thread. Divide by double the number of allocated threads to provide
     * buffer. Samtools' memory limit is wildly inaccurate and unpredictable, in my experience.
     */
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def mem_per_thread = ( avail_mem ).intdiv( task.cpus * 2 )


    """
    sort_and_index_alignment.sh $final_mini_alignment_file $sample_name $output_format $task.cpus $mem_per_thread
    """
}



process SAMTOOLS_MERGE_AND_INDEX_PROC {
     
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/01-REALIGN_SAMPLES", mode: 'copy')

    label 'SORT_AND_INDEX_SAMPLE_PROC'

    input:
        tuple val(sample_name), path(sample_mini_bams)
        val(align_to_ref_tag)
        val(output_format)

    output:
        path '*.sorted.merged.final.{bam,cram}', emit: final_sample_bam

    script:

    /*
     * Calculate the mem per thread. Divide by double the number of allocated threads to provide
     * buffer. Samtools' memory limit is wildly inaccurate and unpredictable, in my experience.
     */
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def mem_per_thread = ( avail_mem ).intdiv( task.cpus * 2 )

    """
    merge_and_index_mini_alignments.sh \$PWD $sample_name $align_to_ref_tag $output_format $task.cpus $mem_per_thread
    """
final_sample_output_file="${sample_name}.${ref_tag}.sorted.merged.final.${out_format}"
   
}



process REALIGN_SAMPLES_PROC {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/01-REALIGN_SAMPLES", mode: 'copy')

	label 'REALIGN_SAMPLES'

	input:
		path(sample_input_file)
		val(align_to_ref)
		val(align_to_ref_tag)
		val(original_ref)
        val(output_format)


	output:
        path '*.{bam,cram}', emit: final_alignment
        path '*.{bam,cram}.*', emit: final_alignment_index

	script:

    def avail_mem = task.memory ? task.memory.toGiga() : 0

    /*
     * Calculate the mem per thread. Divide by double the number of allocated threads to provide
     * buffer. Samtools' memory limit is wildly inaccurate and unpredictable, in my experience.
     */
    def mem_per_thread = ( avail_mem ).intdiv( task.cpus * 2 )

	"""
    echo "Sample file: $sample_input_file"
	bash realign_bwa.sh ${sample_input_file} ${align_to_ref} ${align_to_ref_tag} ${original_ref} ${output_format} $task.cpus $mem_per_thread
	"""
}
