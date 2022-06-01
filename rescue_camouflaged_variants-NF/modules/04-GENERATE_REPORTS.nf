
/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

workflow GENERATE_REPORTS_WF {

    take:
        combined_stats_ch

    main:
	combined_stats_ch.map { println it[0] } 
        combined_stats_ch.map { println it[1] } 

        combined_stats_ch.map { tuple it[0], it[1]  } | view() | set{ sample_tuple } 
        println sample_tuple

	GENERATE_REPORT_PROC( sample_tuple )

}

process GENERATE_REPORT_PROC {

    publishDir("${params.results_dir}/04-REPORT", mode: 'copy')

    label 'GENERATE_REPORT'

    input:
        tuple val(sample_name),
            path('low_depth-merged')
            //path('*.low_mapq-merged.bed'),
            //path('*_perRegionMetrics.bed'),
            //path('*.all.dark.regions.bed'),
            //path('*_extractionCoverage.bed.gz'),
            

    output:
	tuple val(sample_name),
            path("*.html"), emit: report

    script:


    """
	
        runDRFPath=`realpath '${params.results_dir}/01-RUN_DRF'`
        calcStatsPath=`realpath '${params.results_dir}/02-CALCULATE_BAM_STATS'`
        GOI=`realpath '${params.genes_of_interest}'`

	Rscript -e "rmarkdown::render('${params.Report_Rmd}', params=list('runDRFDir'='${runDRFPath}', 'calculateBamStats' = '${calcStatsPath}', 'genesOfInterest' = '${GOI}') , clean=T)"

    """

}



