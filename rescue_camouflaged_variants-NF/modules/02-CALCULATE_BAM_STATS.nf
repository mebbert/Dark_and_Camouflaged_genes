
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

workflow CALCULATE_BAM_STATS_WF {

    take:
	sample_DRF_output_ch

    main:
        
        sample_DRF_output_ch.map { println it[0] } 
        sample_DRF_output_ch.map { println it[1] } 

        sample_DRF_output_ch.map { tuple it[0], it[1]  } | view() | set{ sample_tuple } 
        println sample_tuple

        //sample_DRF_output_file_test = file(sample_DRF_output_file)
//	println sample_DRF_output_file.view()
	//sample_tuple = tuple , sample_DRF_output_file.

	//println sample_tuple
        GENERATE_LOW_MAPQ_AND_DEPTH_BEDS_PROC( sample_tuple )
            | MERGE_DARK_REGIONS_PROC
//            | GENERATE_REPORT_PROC
           // | CALC_BAM_STATS_PROC
	

    emit:
        MERGE_DARK_REGIONS_PROC.out.collect()
}

process GENERATE_LOW_MAPQ_AND_DEPTH_BEDS_PROC {

    publishDir("${params.results_dir}/02-CALCULATE_BAM_STATS", mode: 'copy')

    label 'GENERATE_LOW_MAPQ_AND_DEPTH_BEDS'

    input:
        tuple val(sample_name), val(sample_DRF_output_file)

    output:
        tuple val(low_depth),
                path('*.dark.low_depth.bed.gz'),
                emit: low_depth_out
        tuple val(low_mapq),
                path('*.dark.low_mapq.bed.gz'),
                emit: low_mapq_out

    script:
    
    /*
     * Define low_depth and low_mapq out files
     */
   

    low_depth = sample_DRF_output_file.name.replaceFirst('.bed.gz', '.dark.low_depth.bed.gz')
    low_mapq = sample_DRF_output_file.name.replaceFirst('.bed.gz', '.dark.low_mapq.bed.gz')


    """
    bash combine_DRF.sh "${sample_DRF_output_file}" ${sample_name} $task.cpus
    """

    //combine_DRF_output.py "${sample_DRF_output_file}" "${low_depth}" "${low_mapq}"
}

process MERGE_DARK_REGIONS_PROC {
    
    publishDir("${params.results_dir}/02-CALCULATE_BAM_STATS", mode: 'copy')

    label 'MERGE_DARK_REGIONS'

    input:
        /*
         * This file will be either the file containing low-depth 
         * (i.e., 'dark-by-depth') or low-mapq (i.e., 'dark-my-MAPQ') regions
         */
        tuple val(sample_name), path(low_depth_file)
        tuple val(sample_name), path(low_mapq_file)

    output:
        /*
         * The output file will either be <prefix>.low_depth-merged.bed or
         * prefix.low_mapq-merged.bed
         */
        //tuple val(sample_name),
        //        path('*.low_depth-merged.bed'),
        //        path('*.low_mapq-merged.bed'),
        //        path(dark_region_file),
        //        emit: dark_region_files
        tuple val(sample_name),
                path('*.low_depth-merged.bed'),
                path('*.low_mapq-merged.bed'),
		path('*_perRegionMetrics.bed'),
		path('*.all.dark.regions.bed'),
		path('*_extractionCoverage.bed.gz'),
                emit: dark_region_files

    script:

    low_depth_merged_file = sample_name.replaceFirst("low_mapq.bed.gz", "low_depth-merged.bed")
    low_mapq_merged_file = sample_name.replaceFirst("low_mapq.bed.gz", "low_mapq-merged.bed")
    dark_region_file = sample_name.replaceFirst("low_mapq.bed.gz", "all.dark.regions.bed")
    

    //low_depth_merged_file = low_depth_file.contains("low_depth.bed.gz") ?
    //                dark_region_file.replaceFirst("low_depth.bed.gz", "low_depth-merged.bed") :
    //                dark_region_file.replaceFirst("low_mapq.bed.gz", "low_mapq-merged.bed")
    //low_mapq_merged_file = low_mapq_file.contains("low_mapq.bed.gz") ?
    //                dark_region_file.replaceFirst("low_depth.bed.gz", "low_depth-merged.bed") :
    //                dark_region_file.replaceFirst("low_mapq.bed.gz", "low_mapq-merged.bed")

                //echo "ERROR (`date`): Failed to merge coordinates for ${dark_region_file}. See log for details."
		//echo "ERROR (`date`): Failed to merge coordinates for ${dark_region_file}. See log for details."
    """

    pigz -c -d ${low_mapq_file} | \\
        awk 'NR>1{print \$0}' | remove_unassembled_contigs.py | \\
        bedtools intersect -sorted -wao -a "${params.extraction_bed}" -b - | \\
        pigz -c - > ${sample_name}.low_mapq_extractionCoverage.bed.gz

    pigz -c -d ${low_depth_file} | \\
        awk 'NR>1{print \$0}' | remove_unassembled_contigs.py | \\
        bedtools intersect -sorted -wao -a "${params.extraction_bed}" -b - | \\
        pigz -c - > ${sample_name}.low_depth_extractionCoverage.bed.gz


    pigz -c -d ${low_depth_file} | \\
        awk 'NR>1{print \$0}' | \\
        remove_unassembled_contigs.py | \\
        bedtools intersect -sorted -a "${params.extraction_bed}" -b - -loj | \\
        bedtools merge -i - -o distinct,mean,median,stdev -c 4,10,10,10 > ${sample_name}_low_depth_perRegionMetrics.bed


    pigz -c -d ${low_mapq_file} | \\
        awk 'NR>1{print \$0}' | \\
        remove_unassembled_contigs.py | \\
        bedtools intersect -sorted -a ${params.extraction_bed} -b - -loj | \\
        bedtools merge -i - -o distinct,mean,median,stdev -c 4,10,10,10 > ${sample_name}_low_mapq_perRegionMetrics.bed


    ##################################################################################
    # Merge coordinates for depth bed, removing regions that are less than 20bp long #
    ##################################################################################
    # Swallow 141 (SIGPIPE) because remove_unassembled_contigs.py closes pipe abruptly.
    time bedtools merge -d 20 -c 5 -o mean,median -i ${low_depth_file} | \\
        remove_unassembled_contigs.py | \\
        awk '{ if(\$3 - \$2 > 20) print \$0}' \\
        > ${low_depth_merged_file} \\
        || if [[ \$? -eq 141 ]]; then  # Swallow 141 (SIGPIPE)
                true
            else
                echo "ERROR (`date`): Failed to merge coordinates for ${low_depth_file}. See log for details."
                exit 1
            fi


    time bedtools merge -d 20 -c 5 -o mean,median -i ${low_mapq_file} | \\
        remove_unassembled_contigs.py | \\
        awk '{ if(\$3 - \$2 > 20) print \$0}' \\
        > ${low_mapq_merged_file} \\
        || if [[ \$? -eq 141 ]]; then  # Swallow 141 (SIGPIPE)
                true
            else
                echo "ERROR (`date`): Failed to merge coordinates for ${low_mapq_file}. See log for details."
                exit 1
            fi
  
    time cat ${low_depth_merged_file} ${low_mapq_merged_file} | sort -k1,1 -k2,2n > ${dark_region_file}

    """
}

process GENERATE_REPORT_PROC {

    publishDir("${params.results_dir}/04-REPORT", mode: 'copy')

    label 'GENERATE_REPORT'

    input:
        //tuple val(sample_name),
           // path('low_depth-merged') from dark_region_files.collect()
            //path('*.low_mapq-merged.bed'),
            //path('*_perRegionMetrics.bed'),
            //path('*.all.dark.regions.bed'),
            //path('*_extractionCoverage.bed.gz'),
            
            val info from dark_region_files.collect()

    output:
	tuple val(sample_name),
            path("*.html"), emit: report

    script:


    """
	
        runDRFPath=`realpath '${params.results_dir}/01-RUN_DRF'`
        calcStatsPath=`realpath '${params.results_dir}/02-CALCULATE_BAM_STATS'`
        GOI=`realpath '${params.genes_of_interest}'`
        
        echo \$runDRFPath
        echo \$calcStatsPath
        echo \$GOI

	Rscript -e "rmarkdown::render('${params.Report_Rmd}', params=list('runDRFDir'='\$runDRFPath', 'calculateBamStats' = '\$calcStatsPath', 'genesOfInterest' = '\$GOI') , clean=T, output_file = './${params.masked_ref_tag}_${params.sample_input_tag}_Report.html', output_dir = './')"

    """

}



process CALC_BAM_STATS_PROC {

    publishDir("${params.results_dir}/02-CALCULATE_BAM_STATS/BamStats", mode: 'copy')

    label 'CALC_STATS'

    input:
        tuple val(sample_name),
                 path(merged_low_depth_regions_file),
                 path(merged_low_mapq_regions_file),
                 path(per_region_metric_file),
                 path(full_dark_region_file)
			

    output:
        path("*low_depth_camo_regions.bed")
        path("*low_mapq_camo_regions.bed")
        path("*bam_metrics.txt")


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

    bedtools intersect \\
        -a "${params.extraction_bed}" \\
        -b "${merged_low_mapq_regions_file}" \\
        -loj > "${sample_name}.low_mapq_camo_regions.bed"

    awk -F"\t" 'BEGIN{OFS="\t"}\$6!="."{a[\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$6"\t"\$7"\t"\$8] += \$9; b[\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$6"\t"\$7"\t"\$8] +=1}END{for(x in a){print x,a[x]/b[x]}}' "${sample_name}".low_mapq_camo_regions.bed > "${sample_name}".low_mapq_camo_regions.meanDepth.bed

    awk -F"\t" 'BEGIN{OFS="\t"}\$6!="."{a[\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$6"\t"\$7"\t"\$8] += \$9; b[\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$6"\t"\$7"\t"\$8] +=1}END{for(x in a){print x,a[x]/b[x]}}' "${sample_name}".low_depth_camo_regions.bed > "${sample_name}".low_depth_camo_regions.meanDepth.bed

    calc_bam_metrics.py \\
        "${full_dark_region_file}" \\
        "${merged_low_depth_regions_file}" \\
        "${merged_low_mapq_regions_file}" \\
        > "${sample_name}.bam_metrics.txt"

    """
}




