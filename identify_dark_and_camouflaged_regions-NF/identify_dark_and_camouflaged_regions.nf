// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2

// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##

/*
 * TODO: Write good comments describing these input values.
 *
 * This is the reference the sample(s) is/are *currently* aligned to.
 *
 * NOTE: All input samples must already be aligned to the same reference.
 */
params.original_ref = "${projectDir}/../sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" //add a default value

/*
 * This is the reference genome that samples will be re-aligned to and that will
 * be used to define 'dark' and 'camouflaged' regions.
 */
params.align_to_ref = "${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

/*
 * TODO: Make note that the user *can* specify *.cram or *.bam, but they can
 * also just provide the folder path and we'll take everything with one of those
 * extensions.
 */
// params.input_sample_folder = "${projectDir}/../samples/ADSP/input_sample_folder/*.cram"
// params.input_sample_folder = "${projectDir}/test_data/ADSP_sample_input_sample_folder/*"
params.input_sample_folder = "${projectDir}/original_ADSP_samples/"

/*
 * TODO: Make the pipeline follow this setting?
 */
params.output_format = 'bam'

/*
 *
 */
params.ref_tag = 'hg38'

/*
 *
 */
params.gff = "${projectDir}/../sequencing_resources/annotations/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.93.gff3"

/*
 *
 */
params.sequencer = 'illuminaRL100'

/*
 *
 */
params.results_dir = './results'

/*
 *
 */
params.prefix = 'IlluminaRL100.hg38.combined'

/*
 * Prefix for masked reference files.
 */
params.mask_ref_prefix = 'hg38_camo_mask'

log.info """\
 CALCULATE DARK AND CAMO REGIONS PIPELINE
 ==========================================
 cram reference         : ${params.original_ref}
 reference              : ${params.align_to_ref}
 input_sample_folder    : ${params.input_sample_folder}
 ref_tag                : ${params.ref_tag}
 reference indexes      : ${params.ref_index}
 """

// Import Modules
include {REALIGN_SAMPLES} from './modules/REALIGN_SAMPLES.nf'
include {RUN_DRF} from './modules/RUN_DRF.nf'
include {COMBINE_DRF_OUTPUT} from './modules/COMBINE_DRF_OUTPUT.nf'
//include {step_03_CALC_BAM_METRICS} from './modules/03_CALC_BAM_METRICS.nf'
include {PREPARE_ANNOTATION_BED} from './modules/PREPARE_ANNOTATION_BED.nf'
include {CREATE_BED_FILE} from './modules/CREATE_BED_FILE.nf'

// Define initial files and channels
/*
 * Return only *.sam, *.bam, and *.cram files (no index or other files that may
 * be in the folder).
 *
 * TODO: Make the pipeline actually support .bam files.
*/
samples = Channel.fromPath(params.input_sample_folder)
                    .filter( ~/.*(\.sam|\.bam|\.cram)/ )
align_to_ref = file(params.align_to_ref)
original_ref = file(params.original_ref)
ref_tag = params.ref_tag
run_drf_jar = file("${projectDir}/bin/DarkRegionFinder.jar")
sequencer = params.sequencer


workflow{
	REALIGN_SAMPLES(samples, align_to_ref, params.ref_tag, original_ref)
	//step_00_GET_SAMPLES.out.final_bams.view()
	RUN_DRF(REALIGN_SAMPLES.out.final_bams, align_to_ref, ref_index, run_drf_jar)
	//step_01_RUN_DRF.out.low_mapq_bed.view()
	COMBINE_DRF_OUTPUT(RUN_DRF.out.low_mapq_bed.collect(), params.prefix)
//	step_03_CALC_BAM_METRIC()
	PREPARE_ANNOTATION_BED(params.gff)
	CREATE_BED_FILE(COMBINE_DRF_OUTPUT.out.low_depth_out,
                     COMBINE_DRF_OUTPUT.out.low_mapq_out, align_to_ref, ref_index,
                     PREPARE_ANNOTATION_BED.out.prepped_anno_bed, sequencer,
                     ref_tag)
}
