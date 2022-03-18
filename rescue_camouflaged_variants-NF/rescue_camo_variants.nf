import java.text.SimpleDateFormat

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

/*
 * Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##
 */


/*
 * Get current date and time
 */
def date = new Date()
def sdf = new SimpleDateFormat("yyyy_MM_dd-HH.mm.ss")
def time_stamp = sdf.format(date)

/*
 * The original and UNMASKED reference that the samples are currently aligned to.
 * This will be used to identify false-positive variants in the rescued variant
 * VCF files. These false positives are the result of 'reference-based artifacts',
 * which we describe in the paper and in notes for the respective NextFlow process.
 */
params.unmasked_ref_fasta = "${projectDir}/../references/GRCh38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"

/*
 * The MASKED reference to be *RE-ALIGNED TO*; this should be the same reference
 * as the `unmasked_ref_fasta`, except it needs to have been prepared (i.e., masked)
 * using our pipeline to identify dark and camouflaged regions. If the reference
 * provided here has not been masked for camouflaged regions, we cannot rescue camouflaged
 * variants.
 */
params.masked_ref_fasta = "${projectDir}/results/MASK_GENOME/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fa"

/*
 * A 'tag' to describe the reference genome. This is used for naming output
 * in the NextFlow process that identifies false positives (a.k.a.
 * reference-based artifacts).
 */
params.masked_ref_tag = 'NCBI_GRCh38_no_alt_plus_hs38d1_analysis_set'


/*
 * Path to input bams for rescuing camouflaged variants. Must be an absolute
 * path (not relative).
 */
params.input_sample_path = "${projectDir}/test_data/UKY_ADSP_crams"

/*
 * This string can be used to help name the results folder (if provided and included when defining
 * the results parameter). 
 */
params.sample_input_tag = "Test_samples"


/*
 * Path to the .bed file defining regions to extract reads from for samples
 * (i.e., the camouflaged regions). This MUST contain coordinates based on the
 * reference the samples are *currently* aligned to. It matters where the
 * reference was obtained from (i.e., NCBI, Ensembl, UCSC, etc.); they are NOT
 * interchangeable! The camouflaged regions will change depending on which
 * contigs are included in the reference. There are also other more minor
 * (but super annoying) differences between builds (e.g., chromsome naming
 * and even sorting).
 *
 * We provide .bed files for the following human reference genomes:
 *   1. TODO: List reference genomes we provide camo beds for.
 *
 * If working with a reference genome (from any organism) that we have not
 * already identified camouflaged regions for, it can be prepared using our
 * pipeline for identifying camouflaged regions (de novo).
 */
params.extraction_bed = "${projectDir}/test_data/CR1-extraction-1KG_ref.bed"

/*
 * Path to the 'align_to' .bed file. This is used to focus on just the camoflaged
 * regions in the NextFlow process that identifies false-positive variants.
 * This .bed file comes from the `CREATE_BED_FILE` process in the
 * `Identify Dark and Camouflaged Regions` NextFlow workflow. This .bed file
 * must be specific to the reference genome being used in the `unmasked_ref_fasta`
 * argument.
 */
params.align_to_bed = "${projectDir}/test_data/illuminaRL100.hg38.camo.align_to.sorted.bed"

/*
 * Path to the .bed file that GATK will use to call variants. This MUST
 * contain coordinates based on the *MASKED* reference the samples *WILL BE*
 * aligned to in this pipeline. This can be the same reference *VERSION* the
 * samples were previously aligned to, but it must be masked. This .bed file
 * is similar to the 'align_to.bed' file used to mask the genome, except it is
 * restricted to the coding (CDS) regions within the camouflaged region,
 * rather than the entire camouflaged region.
 *
 * We provide the appropriate GATK .bed files for the same human reference
 * genomes as the extraction_bed.
 *
 * TODO: all camo or just CDS? Was previously only CDS.
 */
params.gatk_bed = "${projectDir}/test_data/CR1-GATK-1KG_ref.bed"

/*
 * 
 */
params.camo_annotations = "${projectDir}/test_data/illuminaRL100.hg38.camo.align_to.sorted.bed"

/*
 * Define the number of samples to run in a single rescue batch. This can
 * be a reasonably large number because each rescue batch is dealing with a
 * relatively small portion of the genome, so the anlysis runs fairly quickly.
 */
params.n_samples_per_batch = 50

/*
 * Path where final results should be saved. The workflow will generate subfolders for each of the
 * steps below with their respective outputs.
 *
 * IMPORTANT: NextFlow is picky about relative paths. i.e., simply using 'results' is not allowed.
 * It must be either an absolute path or a relative path that includes a preceing './' or the global
 * variable "${projectDir}" (e.g., '${projectDir}/').
 */
params.results_dir = "./results/${params.masked_ref_tag}_${params.sample_input_tag}-${time_stamp}"

/*
 * The max number of repeat regions to rescue. Some genomic sequences are
 * repeated many times throughout the genome. We do not currently have much
 * confidence calling variants when there are >5 identical regions; this is
 * really dependent on sequencing depth and is an area for further research.
 */
params.max_repeats_to_rescue = 5

/*
 * Parameter defining whether to clean tmp files throughout the run. Depending
 * on how many samples are being run, this can generate tens-of-thousands of files,
 * though they are generally relatively small.
 */
params.clean_tmp_files = false 

/*
 * Defines the size of intervals to break the align_to_ref into for
 * running/parallelizing DRF.
 *
 * A DRF interval length of 10_000_000 will create nearly 500 jobs *per sample*
 * for the human genome (assuming small contigs are included). Small, unplaced
 * contigs and HLA haplotypes make up nearly 200 of the ~500. This is based on
 * the reference the input samples were originally aligned to ('original_ref').
 */
params.DRF_interval_length = 5_000_000

/*
 * Specify where the DRF jar file can be found
 */
params.DRF_jar = file("${projectDir}/bin/DarkRegionFinder.jar")




log.info """\
 RESCUE CAMO VARIANTS PIPELINE
 ==========================================
 unmasked reference                : ${params.unmasked_ref_fasta}
 masked reference                  : ${params.masked_ref_fasta}
 masked reference tag              : ${params.masked_ref_tag}
 sample input path                 : ${params.input_sample_path}
 sample input tag                  : ${params.sample_input_tag}
 extraction bed                    : ${params.extraction_bed}
 GATK bed                          : ${params.gatk_bed}
 align to bed                      : ${params.align_to_bed}
 camo annotations                  : ${params.camo_annotations}
 n samples per batch               : ${params.n_samples_per_batch}
 results dir                       : ${params.results_dir}
 max repeats to rescue             : ${params.max_repeats_to_rescue}
 clean tmp files?                  : ${params.clean_tmp_files}
 DRF_interval_length               : ${params.DRF_interval_length}
 DRF_jar                           : ${params.DRF_jar}
 """

/*
 * Import Modules
 */
include {RUN_DRF_WF} from './modules/01-RUN_DRF.nf'
include {RESCUE_CAMO_VARS_WF} from './modules/03-RESCUE_CAMO_VARS_PROCS.nf'


workflow{

    RUN_DRF_WF( params.DRF_interval_length )

    CALCULATE_BAM_STATS( RUN_DRF_WF.out )

    RESCUE_CAMO_VARS_WF()

}

