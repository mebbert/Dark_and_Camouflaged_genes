import groovy.io.FileType

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2


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

    publishDir("${params.results_dir}/03-IDENTIFY_FALSE_POSITIVES", mode: 'copy')

	label 'IDENTIFY_FALSE_POSITIVES'

    output:
        path "reference_based_artifacts.${params.masked_ref_tag}_2.bed"

    script:
    """
    artifacts_bed="reference_based_artifacts.${params.masked_ref_tag}_2.bed"

    bash identify_false_positives.sh \\
        "${params.camo_annotations}" \\
        "${params.rescue_gene_elements}" \\
        "${params.mask_bed}" \\
        "${params.unmasked_ref_fasta}" \\
        "\${artifacts_bed}" \\
        "${task.cpus}"
    """
}
