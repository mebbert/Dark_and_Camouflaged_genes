/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

/*
 * Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##
 */


/*
 * TODO: Need to set up proper environment (e.g., load java modules automatically).
 */
 
/*
 * Path to input bams for rescuing camouflaged variants. Must be an absolute
 * path (not relative).
 */
params.bam_path = "${projectDir}/test_data/UKY_ADSP_crams"

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
 * The MASKED reference to be *ALIGNED TO*. If the reference provided here has
 * not been masked for camouflaged regions, we cannot rescue camouflaged
 * variants. This reference must be prepared using our pipeline to identify
 * camouflaged regions.
 */
params.masked_ref_fasta = "${projectDir}/results/MASK_GENOME/Homo_sapiens.GRCh38.dna.primary_assembly-MASKED.fa"

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
 */
params.gatk_bed = "${projectDir}/test_data/CR1-GATK-Ensembl.bed"

/*
 * Define the number of samples to run in a single rescue batch. This can
 * be a reasonably large number because each rescue batch is dealing with a
 * relatively small portion of the genome, so the anlysis runs fairly quickly.
 */
params.n_samples_per_batch = 50

/*
 * Define the directory to publish final results in.
 */
params.results_dir = "./results"

/*
 * The max number of repeat regions to rescue. Some genomic sequences are
 * repeated many times throughout the genome. We do not currently have much
 * confidence calling variants when there are >5 identical regions. This is
 * and area for further research.
 */
params.max_repeats_to_rescue = 5

/*
 * TODO: What was this for?
 */
params.ref_tag = 'hg38'


// TODO: What's this for?
log.info """\
 RESCUE CAMO VARIANTS PIPELINE
 ==========================================
 cram reference         : ${params.cram_ref}
 crams                  : ${params.crams}
 ref_tag                : ${params.ref_tag}
 """

/*
 * Import Modules
 */
include {MASK_GENOME} from './modules/MASK_GENOME.nf'
include {RESCUE_CAMO_VARS_WF} from './modules/RESCUE_CAMO_VARS_PROCS.nf'


/*
 * Define initial files and channels
 */
bam_path = params.bam_path
extraction_bed = file(params.extraction_bed)
masked_ref_fasta = params.masked_ref_fasta
gatk_bed = file(params.gatk_bed)
n_samples_per_batch = params.n_samples_per_batch
max_repeats_to_rescue = params.max_repeats_to_rescue

workflow{
/*
    MASK_GENOME(extraction_bed,
            "${projectDir}/../references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            "${projectDir}/../references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
            'Homo_sapiens.GRCh38.dna.primary_assembly-MASKED')
*/
    RESCUE_CAMO_VARS_WF(bam_path, extraction_bed, masked_ref_fasta, gatk_bed, 
                         n_samples_per_batch, max_repeats_to_rescue)
}

