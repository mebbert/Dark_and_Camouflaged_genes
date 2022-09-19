import groovy.io.FileType

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2


workflow VARIANT_FILTERING_WF {

    take:
        vcf_tuples_ch

    main:

        CAT_VARIANTS_PROC(vcf_tuples_ch)
        COMBINE_PLOIDY_PROC(CAT_VARIANTS_PROC.out.ploidy_vcf.collect())
        GET_VARIANT_METRICS_PROC(COMBINE_PLOIDY_PROC.out.combined_vcf)
        GET_GENE_NUM_PROC(GET_VARIANT_METRICS_PROC.out.variant_metrics)

    emit:
	GET_GENE_NUM_PROC.out

}


/*
 */
process CAT_VARIANTS_PROC {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/01-CAT_VARIANTS", mode: 'copy')

	label 'CAT_VARIANTS'

    input:
        tuple val(ploidy), path(gvcfs)

	output:
		path("*${ploidy}.vcf", emit: ploidy_vcf)

	script:
	def inputMap = [:]
	for(gvcf in gvcfs) {
	    groups = gvcf =~ /full_cohort\.combined\.(\w*)_\d*-\d*.ploidy_\d*\.vcf/
	    chrom = groups[0][1]
	    inputMap[chrom] = (inputMap[chrom] == null) ? "I=$gvcf" : inputMap[chrom] + " I=$gvcf"
	}
	inputFiles = ""
	chromOrder = ["1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y"]
	for(chrom in chromOrder) {
	    inputFiles = (inputMap[chrom] == null) ? inputFiles : "$inputFiles " + inputMap[chrom]
	}

	"""
	    echo $inputMap
	    echo $inputFiles
            out="full_cohort.ADSP.camo_genes.genotyped.hg38.${ploidy}.vcf"
            gatk GatherVcfs\
                O=\$out \
                ${inputFiles} \
                R=${params.unmasked_ref_fasta}
            echo "Completed gatherVCFs"

	"""
}

process COMBINE_PLOIDY_PROC {
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/01-CAT_VARIANTS", mode: 'copy')

        label 'COMBINE_PLOIDY'

    input:
        val(ploidy_vcfs)

        output:
                path('*combined.vcf', emit: combined_vcf)

        script:
        def avail_mem = task.memory ? task.memory.toGiga() : 0

        def Xmx = avail_mem >= 8 ? "-Xmx${avail_mem - 3}G" : ''
        def Xms = avail_mem >= 8 ? "-Xms${avail_mem.intdiv(2)}G" : ''
	def xxFlaf = avail_mem >= 8 ? "-XX:+UseSerialGC" : ''

        """
        input_string="${ploidy_vcfs.join(' I=')}"
        input_string="I= \${input_string}"
	echo ${Xmx}
	echo ${Xms}

        java "${Xmx}" "${Xms}" "${xxFlag}" -jar \$PICARD SortVcf\
            \${input_string} \
            O=full_cohort.ADSP_WES.camo_genes.genotyped.hg38.combined.vcf

        """
}


process GET_VARIANT_METRICS_PROC {
     /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/02-GET_VARIANT_METRICS", mode: 'copy')

        label 'GET_VAR_METRICS'

    input:
        path(combined_vcf)

        output:
                path('filtered.variant_metrics.txt', emit: variant_metrics)

        script:
        """
            ANNOTATED=annotated.vcf
            FILTERED=filtered.vcf
            VARIANT_METRICS=filtered.variant_metrics.txt

            calc_genotype_annotations.py ${combined_vcf} > \$ANNOTATED
            # filter.py \$ANNOTATED ${params.ref_based_artifact_bed} > \$FILTERED
            # extractVariantQualityMetrics.py \$FILTERED > \$VARIANT_METRICS
            extractVariantQualityMetrics.py \$ANNOTATED > \$VARIANT_METRICS

        """
}


process GET_GENE_NUM_PROC {
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/03-GET_GENE_NUM", mode: 'copy')

        label 'GET_GENE_NUM'

    input:
        path(variant_metrics)

        output:
                path('geneNumber*.txt', emit: gene_number_files)

        script:
        """
            awk '\$6 == "PASS" && \$9 >= 2.0' ${variant_metrics} | \
                bedtools intersect \
                    -a ${params.annotation_bed} \
                    -b - | \
                awk '\$4 == "CDS" {print \$NF}' | \
                getGeneNumber.py \
                > geneNumber.QD2.0.txt &

            awk '\$6 == "PASS" && \$9 >= 3.0' ${variant_metrics} | \
                bedtools intersect \
                    -a ${params.annotation_bed} \
                    -b - | \
                awk '\$4 == "CDS" {print \$NF}' | \
                getGeneNumber.py \
                > geneNumber.QD3.0.txt

        """


}



