
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2


/*
 * Run DRF for each individual sample. This workflow is parallelized by splitting
 * the genome into intervals. The results for each sample are then combined to
 * be passed on to the next stage of the workflow which will combine *across*
 * samples.
 */
workflow RUN_DRF_WF {
    take:
        sample_input_ch
        align_to_ref
        DRF_jar
        interval_length

    main:

        /*
         * Create intervals to split DRF jobs across intervals
         */
        // intervals = Channel.from( create_intervals( align_to_ref, interval_length ) )
        intervals_ch = Channel.of('1:10000-20000', '5:55555-66666')



        /*
         * Run DRF on sample for each interval
         */
        // sample_files = Channel.fromPath(sample_input_dir + "*.bam")
        // sample_files = Channel.fromPath(sample_input_dir + "A-CUHS-CU000208-BL-COL-56227BL1.Ensembl_GRCh38_release_93.bam")
        // sample_files.combine(intervals_ch).set { files_and_intervals }
        // RUN_DRF_PROC( files_and_intervals, align_to_ref, DRF_jar)

        /*
         * Create cartesian product for samples and inputs so all intervals
         * are run on all samples
         */
        samples_and_intervals = sample_input_ch
            .combine(intervals_ch)

        RUN_DRF_PROC( samples_and_intervals, align_to_ref, DRF_jar)



        /*
         * Combine sample DRF files
         */

        /* Collect all output from RUN_DRF_PROC and group by sample */
        DRF_output_by_sample = RUN_DRF_PROC.out.low_mapq_beds.collect().flatten()
                        .map {
                            sample_file_path ->
                                // println "sample file path: " + sample_file_path
                                // println "sample file path name: " + sample_file_path.name
                                sample_name = sample_file_path.name.substring(0, sample_file_path.name.indexOf('.'))
                                tuple(
                                        sample_name, sample_file_path
                                     )
                        }
                        .groupTuple()

        COMBINE_SAMPLE_DRF_FILES_PROC(align_to_ref, DRF_output_by_sample)

}



/*
 * Run DRF for a given sample and interval
 */
process RUN_DRF_PROC {
	
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    // publishDir("${params.results_dir}/02-RUN_DRF", mode: 'copy')

	label 'RUN_DRF'

	input:
        tuple path(sample_input_file), val(interval)
//        path(sample_input_file)
//        val(interval)
        val(align_to_ref)
        path(DRF_jar)

	output:
		path '**/*.dark.low_mapq*.bed*.gz', emit: low_mapq_beds
		// path '*_mini_beds/', emit: low_mapq_bed_dir

	script:
	"""
	bash run_DRF.sh ${sample_input_file} ${align_to_ref} ${DRF_jar} ${interval}
	"""
}



/*
 * Combine all of the DRF results across all intervals *within* a given sample.
 */
process COMBINE_SAMPLE_DRF_FILES_PROC {

    publishDir("${params.results_dir}/02-RUN_DRF", mode: 'copy')

    label 'local'

    input:
        val(align_to_ref)

        /* A collection of samples to be staged in current work dir */
        // path(sample_low_mapq_beds)

        tuple val(sample_name), path(sample_low_mapq_beds)

    output:
        path '*final.bed.gz', emit: combined_sample_low_mapq_bed

    script:

    /*
     * Define output file name. Just remove the '.salt_[a-z0-9]+'.
     */
    first_file = sample_low_mapq_beds.first()
    name = first_file.name
    combined_output_file_name = first_file.name.replaceFirst("\\.bed\\.salt_\\w+.gz", ".final.bed.gz")

    
    """
    combine_DRF_sample_output.py $align_to_ref \$PWD $combined_output_file_name
    """
}


/*
 * Define 0-based intervals of size 'interval_length' for every region of the
 * reference genome provided. Any contigs smaller than 'interval_length' will
 * be treated as its own interval.
 */
def create_intervals( ref, interval_length ) {
    
    def ref_faidx = file("${ref}.fai")
    def tmp_interval, remaining_contig_length, intervals = []
    ref_faidx.withReader { reader ->
        while ((line = reader.readLine()) != null) {
            // println "#################"
            // println "${line}"
            toks = line.split("\t")
            contig = toks[0]
            contig_length = toks[1].toInteger()

            remaining_contig_length = contig_length

            /*
             * Split contig into intervals of size 'interval_length'
             */
            // println "Intervals:"
            while (remaining_contig_length > 0) {
                
                /*
                 * 0-based start & end
                 */
                start = 0 + (contig_length - remaining_contig_length)

                if( remaining_contig_length < interval_length ) {
                    end = contig_length
                }
                else {
                    end = start + interval_length
                }
                tmp_interval = "${contig}:${start}-${end}"
                // println tmp_interval
                intervals.add(tmp_interval)

                remaining_contig_length -= interval_length
            }
            // println ""
        }
    }

    return intervals

}
