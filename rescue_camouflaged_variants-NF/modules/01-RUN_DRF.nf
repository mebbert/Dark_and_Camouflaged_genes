
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

    main:

        /*
         * Create sample input channel
         *
         * Recursively collect all .(cr|b)am files (paired with respective index
         * files) from the user-provided input_sample_path.
         *
	 * This should create a tuple of samplename to path for each file passed to this workflow
	 *
	 * TODO: We need to make sure to check that the Crai is available for the cram. If not, we could kick off samtools index or tell the user to do it.
	 * TODO: We will also need to add in some logic to make sure that there is a "/" at the end of the input_sample_path. otherwise it seemingly pulls files at random. 
         */
        Channel.fromFilePairs( "${params.input_sample_path}**.{bam,cram}{,.bai,.crai}", checkIfExists: true )
            | map { tuple(file(it[1][0]).baseName, it[1][0], it[1][1]) }
            | view()
            | set { sample_input_ch }

        /*
         * Create intervals to split DRF jobs across intervals
         */
        intervals = create_intervals( params.current_ref_fasta, params.DRF_interval_length )

	//println intervals.size()
//        intervals_file = file( "${params.current_ref_fasta_tag}.DRF_intervals.txt" )
//        intervals_file.delete()
//        intervals.each {
//            intervals_file << "${it}\n"
//        }

        // intervals = ['1:10000-20000', '1:207496157-207641765',
        //                 '5:55555-66666', '22:15693544-15720708']
        intervals_ch = Channel.from( intervals )


        /*
         * Create cartesian product for samples and inputs so all intervals
         * are run on all samples
         */
        samples_and_intervals = sample_input_ch
            .combine(intervals_ch)
        /*
         * Run DRF and combine sample DRF files
         */
        RUN_DRF_PROC( samples_and_intervals )
            | groupTuple(size: intervals.size())
            | COMBINE_SAMPLE_DRF_FILES_PROC
            | CALC_GLOBAL_STATS_PROC

	CALC_GENES_OF_INTEREST_PROC( sample_input_ch )

    emit:
        COMBINE_SAMPLE_DRF_FILES_PROC.out

}



/*
 * Run DRF for a given sample and interval. The 'sample_input_file' and
 * 'sample_input_file_index' are a (b|cr)am file and its index.
 */
process RUN_DRF_PROC {



    tag { "${sample_name}:${sample_input_file}" }
	
	label 'RUN_DRF_PROC'

	time '1h'

	input:
        tuple val(sample_name),
                path(sample_input_file),
                path(sample_input_file_index),
                val(interval)

	output:
		tuple val(sample_name),
                 path("**/${sample_name}*.dark.low_mapq*.bed*.gz"),
                 emit: low_mapq_beds

	script:
	"""
	bash run_DRF.sh \\
        ${sample_name} \\
        ${sample_input_file} \\
        ${params.current_ref_fasta} \\
        ${params.DRF_jar} \\
        ${interval}
	"""
}



/*
 * Combine all of the DRF results across all intervals *within* a given sample.
 */
process COMBINE_SAMPLE_DRF_FILES_PROC {

    tag { "${sample_name}" }

    publishDir("${params.results_dir}/01-RUN_DRF", mode: 'copy')

    label 'local'

    input:
        /* A collection of samples to be staged in current work dir */
        tuple val(sample_name), path(sample_low_mapq_beds)

    output:
        tuple val(combined_output_file_name),
         path('*final.bed.gz'), emit: combined_low_mapq_bed_ch

//        Channel.fromPath(path('*final.bed.gz'))
//            | map{file(it).baseName, it}
//            | view() 
//            | set{combined_low_mapq_bed_ch}, emit: combined_low_mapq_bed_ch

    script:

    /*
     * Define output file name. Just remove the '.salt_[a-z0-9]+'.
     */
    first_file = sample_low_mapq_beds.first()
    name = first_file.name
    combined_output_file_name = first_file.name.replaceFirst("\\.bed\\.salt_\\w+.gz", ".final.bed.gz")

    
    """
    combine_DRF_sample_output.py \\
        $params.current_ref_fasta \\
        \$PWD \\
        $combined_output_file_name
    """
}

process CALC_GLOBAL_STATS_PROC {

    tag { "${sample_name}" }

    publishDir("${params.results_dir}/01-RUN_DRF", mode: 'copy')

    label 'GLOBAL_STATS_PROC'

    input:
        tuple val(combined_output_file_name),
         path(final_bed)

    output:
       tuple val(combined_output_file_name),
        path('*.GlobalStats.txt'), emit: global_stats

    script:

    combined_output_file_name = combined_output_file_name.replaceFirst("\\.bed.gz", ".GlobalStats.txt")

    """

    bedtools merge -sorted -i ${final_bed} -c 5 -o sum > temp.txt
    calc_perBaseStats.py temp.txt ${combined_output_file_name}
    rm temp.txt

    """
}

process CALC_GENES_OF_INTEREST_PROC {

        tag { "${sample[0]}" }

        publishDir("${params.results_dir}/01-RUN_DRF", mode: 'copy')
	
	label 'CALC_GENES_OF_INTEREST_PROC'


	input:
	 	val(sample)

	output:
		tuple val(sample),
                 path("*_cov.perBase.GOI.bed"),
                 emit: GOI_Coverage

	script:
	"""
 	samtools view -L "${params.genes_of_interest}" -O bam -o "${sample[0]}_ROI.bam" "${sample[1]}"
 	bedtools genomecov -bg -ibam "${sample[0]}_ROI.bam" | \\
 	bedtools intersect -wao -sorted -a - -b "${params.genes_of_interest}" > "${sample[0]}_cov.perBase.GOI.bed"
	

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
                 * 1-based start & end
                 */
                start = 1 + (contig_length - remaining_contig_length)

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
