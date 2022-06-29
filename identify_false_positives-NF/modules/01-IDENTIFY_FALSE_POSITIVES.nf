import groovy.io.FileType

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

workflow IDENTIFY_FALSE_POSITIVES_WF {
    main: 
        IDENTIFY_FALSE_POSITIVES_PROC()
        SPLIT_BLAT_BED_PROC(IDENTIFY_FALSE_POSITIVES_PROC.out.blat_bed)
        EXTRACT_FALSE_POSITIVES_PROC(SPLIT_BLAT_BED_PROC.out.sep_blat_beds.flatMap())
        COMBINE_FALSE_POSITIVE_REGIONS_PROC(EXTRACT_FALSE_POSITIVES_PROC.out.false_positives.collect())
    
    emit:
        COMBINE_FALSE_POSITIVE_REGIONS_PROC.out

}


/*
 * The process to identify false-positive variant calls in the rescued camo
 * variants. These false-positives specifically originate from what we term
 * "reference-based artifacts". These artifacts show up when when two regions
 * in a given camo set are not 100% identical. Thus, when a read from one
 * region is forced to align to the other, a false variant is called.
 *
 * This process creates a list of all possible reference-based artifacts by
 * taking all camo CDS regions with repeat number <= 5 (i.e., <= 5 camo regions
 * in the set) and BLAT them against the whole genome. Any DNA sequence from a
 * hit with sequence identity >= 98% are locally realigned back to the query
 * seqeunce using Bio.pairwise2. Any mismatches or gaps in the aligned sequence
 * are converted into variant positions and output to a bed file that lists
 * positions of reference-based-artifacts to be filtered from the VCF files.
 *
 * NOTE: This only needs to be run once per reference genome.
 */
process IDENTIFY_FALSE_POSITIVES_PROC {

    publishDir("${params.results_dir}/01-IDENTIFY_FALSE_POSITIVES", mode: 'copy')

	label 'IDENTIFY_FALSE_POSITIVES'

    output:
        path "./blat_result/blat.results.bed", emit: blat_bed

    script:
    """

    bash identify_false_positives.sh \\
        "${params.camo_annotations}" \\
        "${params.rescue_gene_elements}" \\
        "${params.mask_bed}" \\
        "${params.unmasked_ref_fasta}" \\
        "${task.cpus}"

    """
}

process SPLIT_BLAT_BED_PROC {
    label 'local'
    input:
        path(blat_bed)
    output:
        path "blat.results.bed.*", emit: sep_blat_beds

    script:
    """
    echo "`date` spliting up the blat_bed"
    
    total_lines=\$(wc -l < ${blat_bed})
    ((lines_per_file = (total_lines + 1000 - 1) / 1000))
    split --lines=\${lines_per_file} ${blat_bed} blat.results.bed.

    """

}

process EXTRACT_FALSE_POSITIVES_PROC {
     
     publishDir("${params.results_dir}/01-IDENTIFY_FALSE_POSITIVES", mode: 'copy')
 
         label 'IDENTIFY_FALSE_POSITIVES'
     input:
        path(sep_blat_beds)
 
     output: 
        path "${sep_blat_beds}.false_positives.txt", emit: false_positives
         
     script:
     """
     echo "`date` extracting false positives"
 
     artifacts_bed="reference_based_artifacts.${params.masked_ref_tag}.bed"
     # create a bed file of positions of the reference base artifacts that need to be filtered from the VCF files
     if ! extract_false_positives.py \\
         "${sep_blat_beds}" \\
         "${params.unmasked_ref_fasta}"  \\
         > "${sep_blat_beds}.false_positives.txt"; then 
         echo "`date` extract_false_positives.py failed for ${sep_blat_beds}"
         exit 1
     fi
     """
}

process COMBINE_FALSE_POSITIVE_REGIONS_PROC {
    publishDir("${params.results_dir}/01-IDENTIFY_FALSE_POSITIVES", mode: 'copy')

    label 'IDENTIFY_FALSE_POSITIVES'
    input:
        path(false_positives)
    output:
        path "reference_based_artifacts.${params.masked_ref_tag}.bed"

    script:
    """
    touch combined_false_positives.txt
    chmod 775 combined_false_positives.txt
    for file in ${false_positives}; do
        cat \$file >> combined_false_positives.txt
    done
    
    bedtools sort -i "combined_false_positives.txt" -g "${params.unmasked_ref_fasta}.fai" > "false_positives_sorted.txt"

    artifacts_bed="reference_based_artifacts.${params.masked_ref_tag}.bed"
    cp "false_positives_sorted.txt" \${artifacts_bed}
    
    """ 

}
