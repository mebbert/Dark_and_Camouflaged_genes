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

}


/*
 */
process CAT_VARIANTS_PROC {

    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/01-CAT_VARIANTS", mode: 'copy')

	label 'local'

    input:
        tuple val(ploidy), path(gvcfs)

	output:
		path('*${ploidy}.vcf', emit: ploidy_vcf)

	script:
	"""
            inputFiles=""
            for gvcf in ${gvcfs}
            do
                echo \${gvcf}
                inputFiles="\${inputFiles} I=\${gvcf}"
            done
            echo \$inputFiles
            out="full_cohort.ADSP.camo_genes.genotyped.hg38.${ploidy}.vcf"
            gatk GatherVcfs\
                O=\$out \
                \$inputFiles \
                R=${params.unmasked_ref_fasta}

	"""
}

process COMBINE_PLOIDY_PROC {
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/01-CAT_VARIANTS", mode: 'copy')

        label 'local'

    input:
        val(ploidy_vcfs)

        output:
                path('*combined.vcf', emit: combined_vcf)

        script:
        """
        input_string=""
        for ploidy_vcf in ${ploidy_vcfs}
        do
            \$input_string="\${input_string} I=\${ploidy_vcf}"
        done

        \$PICARD SortVcf\
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

        label 'local'

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
            filter.py \$ANNOTATED ${params.ref_based_artifact_bed} > \$FILTERED
            extractVariantQualityMetrics.py \$FILTERED > \$VARIANT_METRICS

        """
}


process GET_GENE_NUM_PROC {
    /*
     * Publish results. 'mode: copy' will copy the files into the publishDir
     * rather than only making links.
     */
    publishDir("${params.results_dir}/04-VARIANT_FILTERING/03-GET_GENE_NUM", mode: 'copy')

        label 'local'

    input:
        path(variant_metrics)

        output:
                path('geneNumber*.txt', emit: gene_number_files)

        script:
        """
            awk '\$6 == "PASS" && \$9 >= 2.0' ${variant_metrics} | \
                bedtools interesct \
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



