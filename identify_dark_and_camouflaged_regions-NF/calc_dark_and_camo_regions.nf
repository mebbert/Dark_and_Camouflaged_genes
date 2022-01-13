// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2

// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##

params.cram_ref = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/GRCh38_full_analysis_set_plus_decoy_hla.fa' //add a default value
params.ref = '/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.ref_index='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
params.ref_index_amb='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.amb'
params.ref_index_pac='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.pac'
params.ref_index_ann='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann'
params.ref_index_sa='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.sa'
params.ref_index_bwt='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt'
params.crams = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/ADSP_sample_CRAMs/*.cram'
//params.crams = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/nextflow/test_data/ADSP_sample_crams/*'
params.ref_tag = 'hg38'
params.gff = '/project/mteb223_uksr/sequencing_resources/annotations/hg38_release_93/Homo_sapiens.GRCh38.93.gff3'
params.threads = 16
params.sequencer = 'illuminaRL100'
params.output_dir = './results_dir'
params.prefix = 'IlluminaRL100.hg38.combined'

/*
 * Prefix for masked reference files.
 */
params.mask_ref_prefix = 'hg38_camo_mask'

log.info """\
 CALCULATE DARK AND CAMO REGIONS PIPELINE
 ==========================================
 cram reference         : ${params.cram_ref}
 reference              : ${params.ref}
 crams                  : ${params.crams}
 ref_tag                : ${params.ref_tag}
 reference indexes      : ${params.ref_index}
 """

// Import Modules
include {step_00_GET_BAMS} from './modules/00_GET_BAMS.nf'
include {step_01_RUN_DRF} from './modules/01_RUN_DRF.nf'
include {step_02_COMBINE_DRF_OUTPUT} from './modules/02_COMBINE_DRF_OUTPUT.nf'
//include {step_03_CALC_BAM_METRICS} from './modules/03_CALC_BAM_METRICS.nf'
include {step_04_PREPARE_ANNOTATION_BED} from './modules/04_PREPARE_ANNOTATION_BED.nf'
include {step_05_CREATE_BED_FILE} from './modules/05_CREATE_BED_FILE.nf'

// Define initial files and channels
adsp_samples = Channel.fromPath(params.crams)
ref_index = file(params.ref_index)
ref_index_amb = file(params.ref_index_amb)
ref_index_pac = file(params.ref_index_pac)
ref_index_ann = file(params.ref_index_ann)
ref_index_sa = file(params.ref_index_sa)
ref_index_bwt = file(params.ref_index_bwt)
ref = file(params.ref)
cram_ref = file(params.cram_ref)
ref_tag = params.ref_tag
run_drf_jar = file('/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/Dark_and_Camouflaged_genes/scripts/01_RUN_DRF/DarkRegionFinder.jar')
threads = params.threads
sequencer = params.sequencer
results_dir = params.output_dir


workflow{
	step_00_GET_BAMS(adsp_samples, ref, params.ref_tag, cram_ref, ref_index, ref_index_amb, ref_index_pac, ref_index_ann, ref_index_sa, ref_index_bwt)
	//step_00_GET_BAMS.out.final_bams.view()
	step_01_RUN_DRF(step_00_GET_BAMS.out.final_bams, ref, ref_index, run_drf_jar)
	//step_01_RUN_DRF.out.low_mapq_bed.view()
	step_02_COMBINE_DRF_OUTPUT(step_01_RUN_DRF.out.low_mapq_bed.collect(), params.prefix)
//	step_03_CALC_BAM_METRIC()
	step_04_PREPARE_ANNOTATION_BED(params.gff)
	step_05_CREATE_BED_FILE(step_02_COMBINE_DRF_OUTPUT.out.low_depth_out, step_02_COMBINE_DRF_OUTPUT.out.low_mapq_out, ref, ref_index, step_04_PREPARE_ANNOTATION_BED.out.prepped_anno_bed, sequencer, ref_tag, threads, results_dir)
}
