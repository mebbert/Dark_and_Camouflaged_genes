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
 * The MASKED reference to be *RE-ALIGNED TO*; this could be any reference,
 * except it must have been prepared (i.e., masked) using our pipeline to
 * identify dark and camouflaged regions. If the reference provided here has
 * not been masked for camouflaged regions, rescuing camouflaged variants will
 * not work.
 */
params.masked_ref_fasta = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_04_08-17.34.55/06-MASK_GENOME/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.fa"


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
params.mask_bed = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_04_08-17.34.55/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo.mask_bed.sorted.bed"

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
/*
 * This is the *.camo_annotations.txt file that comes from the step
 * 05-CREATE_BED_FILE in the workflow to define camouflaged regions.
 */
params.camo_annotations = "${projectDir}/../identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_04_08-17.34.55/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo_annotations.txt"


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



log.info """\
 RESCUE CAMO VARIANTS PIPELINE
 ==========================================
 unmasked reference                : ${params.unmasked_ref_fasta}
 masked reference                  : ${params.masked_ref_fasta}
 masked reference tag              : ${params.masked_ref_tag}
 align to bed                      : ${params.mask_bed}
 camo annotations                  : ${params.camo_annotations}
 results dir                       : ${params.results_dir}
 """

/*
 * Import Modules
 */
include {IDENTIFY_FALSE_POSITIVES_PROC} from './modules/01-IDENTIFY_FALSE_POSITIVES.nf'


workflow{

    IDENTIFY_FALSE_POSITIVES_PROC()

}

