
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

workflow CALCULATE_BAM_STATS_WF {

    take:
        sample_DRF_output_file

    main:
        
        sample_DRF_output_file.map { println "############it: $it" }
//        sample_DRF_output_file_ch = Channel.of( file(sample_DRF_output_file) )
//            | map { tuple( it.baseName, it ) }

//        GENERATE_LOW_MAPQ_AND_DEPTH_BEDS_PROC( sample_tuple )
//            | MERGE_DARK_REGIONS_PROC
//            | CALC_BAM_STATS_PROC

    emit:
        MERGE_DARK_REGIONS_PROC.collect()
}

process GENERATE_LOW_MAPQ_AND_DEPTH_BEDS_PROC {

    label 'GENERATE_LOW_MAPQ_AND_DEPTH_BEDS'

    input:
        tuple val(sample_name), path(sample_DRF_output_file)

    output:
        tuple val(sample_name),
                path('*.dark.low_depth.bed.gz'),
                emit: low_depth_out
        tuple val(sample_name),
                path('*.dark.low_mapq.bed.gz'),
                emit: low_mapq_out

    script:
    
    /*
     * Define low_depth and low_mapq out files
     */
    low_depth = sample_DRF_output_file.name.replaceFirst('.bed.gz', 'low_depth.bed.gz')
    low_mapq = sample_DRF_output_file.name.replaceFirst('.bed.gz', 'low_mapq.bed.gz')

    """
    combine_DRF_output.py "${sample_DRF_output_file}" "${low_depth}" "${low_mapq}"
    """
}

process MERGE_DARK_REGIONS_PROC {
    
    label 'MERGE_DARK_REGIONS'

    input:
        /*
         * This file will be either the file containing low-depth 
         * (i.e., 'dark-by-depth') or low-mapq (i.e., 'dark-my-MAPQ') regions
         */
        tuple val(sample_name), path(dark_region_file)

    output:
        /*
         * The output file will either be <prefix>.low_depth-merged.bed or
         * prefix.low_mapq-merged.bed
         */
        tuple val(sample_name),
                path('*.low_depth-merged.bed'),
                path('*.low_mapq-merged.bed'),
                path(dark_region_file),
                emit: dark_region_files

    script:

    merged_file = dark_region_file.contains("low_depth.bed.gz") ?
                    dark_region_file.replaceFirst("low_depth.bed.gz", "low_depth-merged.bed") :
                    dark_region_file.replaceFirst("low_mapq.bed.gz", "low_mapq-merged.bed")

    """

    ##################################################################################
    # Merge coordinates for depth bed, removing regions that are less than 20bp long #
    ##################################################################################
    # Swallow 141 (SIGPIPE) because remove_unassembled_contigs.py closes pipe abruptly.
    time bedtools merge -d 20 -c 5 -o mean,median -i ${dark_region_file} | \\
        remove_unassembled_contigs.py | \\
        awk '{ if(\$3 - \$2 > 20) print \$0}' \\
        > ${merged_file} \\
        || if [[ \$? -eq 141 ]]; then  # Swallow 141 (SIGPIPE)
                true
            else
                echo "ERROR (`date`): Failed to merge coordinates for ${dark_region_file}. See log for details."
                exit 1
            fi

    """
}

process CALC_BAM_STATS_PROC {

    label 'CALC_STATS'

    input:
        tuple val(sample_name),
                 path(merged_low_depth_regions_file),
                 path(merged_low_mapq_regions_file),
                 path(full_dark_region_file)

    output:
        path()

    script:

    """

    # Run bedtools intersect between the low-depth .bed file from this sample
    # and the 'extraction' .bed file of all camouflaged regions identified in 
    # the reference genome version this sample was aligned to. and critical genes .bed
    # to determine if any of the critical genes have dark regions. If so,
    # there may be regions in those genes where variants cannot be reliably
    # rescued for this specific sample (i.e., could be degraded DNA, or intronic
    # region from an exome, etc.).

    bedtools intersect \\
        -a "${params.extraction_bed}" \\
        -b "${merged_low_depth_regions_file}" \\
        -loj > "${sample_name}.low_depth_camo_regions.bed"


    calc_bam_metrics.py \\
        "${full_dark_region_file}" \\
        "${merged_low_depth_regions_file}" \\
        "${merged_low_mapq_regions_file}" \\
        > "${sample_name}.bam_metrics.txt"

    
    """
}
