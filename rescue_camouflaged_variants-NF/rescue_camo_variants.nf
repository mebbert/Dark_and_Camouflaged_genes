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
 * Path to the .fasta the input samples are *currently* aligned to. This is used
 * to run DRF to assess the sample's quality and nature (e.g., if it's exome,
 * whole genome, etc.).
 */
params.current_ref_fasta = "${projectDir}/../references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa"

/*
 * Path to the .bed file defining regions to extract reads from
 * (i.e., the camouflaged regions). This MUST contain coordinates based on the
 * reference the input samples are *currently* aligned to. It matters where the
 * reference was obtained from (i.e., NCBI, Ensembl, UCSC, etc.); they are NOT
 * necessarily interchangeable! While coordinates for a given gene will remain
 * the same for a given genome build (e.g., GRCh38), WHICH regions are
 * camouflaged will change depending on which contigs are included in the
 * reference. There are also other more minor (but super annoying) differences
 * between builds (e.g., chromsome naming and even sorting).
 *
 * We provide .bed files for the following human reference genomes:
 *   1. TODO: List reference genomes we provide camo beds for.
 *
 * If working with a reference genome (from any organism) that we have not
 * already identified camouflaged regions for, it can be prepared using our
 * pipeline for identifying camouflaged regions (de novo).
 */
params.extraction_bed = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_05_17-11.58.12/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo.extraction.sorted.bed"


/*
 * The MASKED reference to be *RE-ALIGNED TO*; this could be any reference,
 * except it must have been prepared (i.e., masked) using our pipeline to
 * identify dark and camouflaged regions. If the reference provided here has
 * not been masked for camouflaged regions, rescuing camouflaged variants will
 * not work.
 */
params.masked_ref_fasta = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_05_17-11.58.12/06-MASK_GENOME/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.fa"


/*
 * The UNMASKED version of the above reference that the samples will be aligned
 * to. This will be used to identify false-positive variants in the rescued variant
 * VCF files. These false positives are the result of what we termed 
 * 'reference-based artifacts', which we describe in the paper and in notes
 * for the respective NextFlow process.
 */
params.unmasked_ref_fasta = "${projectDir}/../references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa"

/*
 * A 'tag' to describe the reference genome. This is used for naming output
 * in the NextFlow process that identifies false positives (a.k.a.
 * reference-based artifacts).
 */
params.masked_ref_tag = 'illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla'


/*
 * Path to input bams for rescuing camouflaged variants. Must be an absolute
 * path (not relative). This path should be to the parent folder for *all*
 * .(cr|b)am files to be included in this rescue run. The workflow will 
 * recursively find and include all .(cr|b)am files in the provided path.
 */
params.input_sample_path = "${projectDir}/test_data/UKY_ADSP_test_crams"

/*
 * This string can be used to help name the results folder (if provided and included when defining
 * the results parameter). 
 */
params.sample_input_tag = "Test_samples"

/*
 * Path to the 'mask.bed' file (previously knwown as the 'align_to' .bed file).
 * This .bed file was used to mask the reference genome, leaving only one unique
 * camouflaged region per set of camouflaged regions. The one unmasked
 * camouflaged remaining has also been expanded by 50bp on each side to allow
 * reads to align at the end of the camouflaged region. 
 *
 * In this pipeline, however, the 'mask.bed' is used to help identify
 * reference-based artifacts (i.e., false-positive variants); it's needed for
 * BLAT. 
 *
 * This .bed file comes from the `05-CREATE_BED_FILE` process in the
 * `Identify Dark and Camouflaged Regions` NextFlow workflow. This .bed file
 * must have coordinates specific to the reference genome being used in the
 * 'masked_ref_fasta' and the `unmasked_ref_fasta` arguments.
 */
params.mask_bed = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_05_17-11.58.12/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo.mask_bed.sorted.bed"

/*
 * Path to the .bed file that GATK will use to call variants. This MUST
 * contain coordinates based on the *MASKED* reference the samples *WILL BE*
 * aligned to in this pipeline. This can be the same reference *VERSION* the
 * samples were previously aligned to, but it must be masked. This .bed file
 * is similar to the 'mask_bed.bed' file used to mask the genome, except it is
 * restricted to the exact camouflaged regions rather than the expanded
 * camouflaged region needed for alignment (i.e., so reads can align to the
 * end of a camouflaged region).
 *
 * This .bed file should be either the *CDS_regions_only* or the *all_camo_regions*
 * .bed file from step 05-CREATE_BED_FILE in the workflow to define camouflaged
 * regions.
 *
 * We provide the appropriate GATK .bed files for the same human reference
 * genomes as the extraction_bed.
 *
 * TODO: Was previously only CDS. Currently adding support for all gene-body
 * elements, but we still ignore anything outside of gene bodies.
 */
params.gatk_bed = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_05_17-11.58.12/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo.GATK.all_camo_regions.bed"

/*
 * This is the *.annotation.bed file that comes from the step
 * 04-PREPARE_ANNOTATION_BED in the workflow to define camo regions. 
 */
params.annotation_bed = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_05_17-11.58.12/04-PREPARE_ANNOTATION_BED/Homo_sapiens.GRCh38.105.chr_added.annotation.bed"

/*
 * This is the *.camo_annotations.txt file that comes from the step
 * 05-CREATE_BED_FILE in the workflow to define camouflaged regions.
 */
params.camo_annotations = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_05_17-11.58.12/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo_annotations.txt"

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
 * Define which gene region elements to rescue variants for. Possible options
 * include "all" or "CDS", where "all" will rescue variants in any gene body
 * element (e.g., including introns and UTR regions) and "CDS will only rescue
 * variants in protein-coding regions.
 *
 * TODO: Make the pipeline also rescue variants *outside* of known gene bodies.
 */
params.rescue_gene_elements = "all"

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
 * This is a bed file of the genes of interest to you. 
 * It's formatted chr\tstart of CDS regions\tend of CDS region\tGene Name 
 */
params.genes_of_interest = "genesOfInterest.bed"

/*
 * Defines the size of intervals to break the reference into for
 * running/parallelizing DRF. This is for the reference that the input samples
 * are currently aligned to.
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

/*
 * Specify where the Report Rmarkdown file is located
 */
params.Report_Rmd = file("${projectDir}/bin/SampleGenomeStats.Rmd")



params.ref_based_artifact_bed = file('/scratch/mlpa241/41/748c8a52ae87f925b06d12fee51b09/reference_based_artifacts.illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.bed')


log.info """\
 RESCUE CAMO VARIANTS PIPELINE
 ==========================================
 current reference                 : ${params.current_ref_fasta}
 unmasked reference                : ${params.unmasked_ref_fasta}
 masked reference                  : ${params.masked_ref_fasta}
 masked reference tag              : ${params.masked_ref_tag}
 sample input path                 : ${params.input_sample_path}
 sample input tag                  : ${params.sample_input_tag}
 extraction bed                    : ${params.extraction_bed}
 GATK bed                          : ${params.gatk_bed}
 align to bed                      : ${params.mask_bed}
 camo annotations                  : ${params.camo_annotations}
 n samples per batch               : ${params.n_samples_per_batch}
 results dir                       : ${params.results_dir}
 max repeats to rescue             : ${params.max_repeats_to_rescue}
 clean tmp files?                  : ${params.clean_tmp_files}
 DRF_interval_length               : ${params.DRF_interval_length}
 Genes_of_interest                 : ${params.genes_of_interest}
 Report_Rmd_Path                   : ${params.Report_Rmd}
 DRF_jar                           : ${params.DRF_jar}
 reference based artifacts         : ${params.ref_based_artifact_bed}
 """

/*
 * Import Modules
 */
include {RUN_DRF_WF} from './modules/01-RUN_DRF.nf'
include {CALCULATE_BAM_STATS_WF} from './modules/02-CALCULATE_BAM_STATS.nf'
include {RESCUE_CAMO_VARS_WF} from './modules/03-RESCUE_CAMO_VARS_PROCS.nf'
include {GENERATE_REPORTS_WF} from './modules/04-GENERATE_REPORTS.nf'
include {VARIANT_FILTERING_WF} from './modules/04-VARIANT_FILTERING.nf'


workflow{

    //RUN_DRF_WF()

    //println RUN_DRF_WF.out
    //CALCULATE_BAM_STATS_WF( RUN_DRF_WF.out ) 

    //println CALCULATE_BAM_STATS_WF.out

    //GENERATE_REPORTS_WF( CALCULATE_BAM_STATS_WF.out  )

    RESCUE_CAMO_VARS_WF()
    VARIANT_FILTERING_WF(RESCUE_CAMO_VARS_WF.out)


}

