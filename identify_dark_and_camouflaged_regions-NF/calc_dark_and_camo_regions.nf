// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2

// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##

params.cram_ref = "${projectDir}/../sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" //add a default value
params.ref = "${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.ref_index="${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
params.ref_index_amb="${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.amb"
params.ref_index_pac="${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.pac"
params.ref_index_ann="${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann"
params.ref_index_sa="${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.sa"
params.ref_index_bwt="${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt"
//params.ref = "${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa"
//params.ref_index="${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa.fai"
//params.ref_index_amb="${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa.amb"
//params.ref_index_pac="${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa.pac"
//params.ref_index_ann="${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa.ann"
//params.ref_index_sa="${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa.sa"
//params.ref_index_bwt="${projectDir}/../sequencing_resources/references/Homo_sapiens.GRCh38_onlyChr22.fa.bwt"

// params.crams = "${projectDir}/../samples/ADSP/crams/*.cram"
params.crams = "${projectDir}/test_data/ADSP_sample_crams/*"
params.ref_tag = 'hg38'
params.gff = "${projectDir}/../sequencing_resources/annotations/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.93.gff3"
params.threads = 16
params.sequencer = 'illuminaRL100'
params.results_dir = './results'
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
include {REALIGN_BAMS} from './modules/REALIGN_BAMS.nf'
include {RUN_DRF} from './modules/RUN_DRF.nf'
include {COMBINE_DRF_OUTPUT} from './modules/COMBINE_DRF_OUTPUT.nf'
//include {step_03_CALC_BAM_METRICS} from './modules/03_CALC_BAM_METRICS.nf'
include {PREPARE_ANNOTATION_BED} from './modules/PREPARE_ANNOTATION_BED.nf'
include {CREATE_BED_FILE} from './modules/CREATE_BED_FILE.nf'

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
run_drf_jar = file("${projectDir}/bin/DarkRegionFinder.jar")
threads = params.threads
sequencer = params.sequencer
results_dir = params.results_dir


workflow{
	REALIGN_BAMS(adsp_samples, ref, params.ref_tag, cram_ref, ref_index, ref_index_amb, ref_index_pac, ref_index_ann, ref_index_sa, ref_index_bwt)
	//step_00_GET_BAMS.out.final_bams.view()
	RUN_DRF(REALIGN_BAMS.out.final_bams, ref, ref_index, run_drf_jar)
	//step_01_RUN_DRF.out.low_mapq_bed.view()
	COMBINE_DRF_OUTPUT(RUN_DRF.out.low_mapq_bed.collect(), params.prefix)
//	step_03_CALC_BAM_METRIC()
	PREPARE_ANNOTATION_BED(params.gff)
	CREATE_BED_FILE(COMBINE_DRF_OUTPUT.out.low_depth_out,
                     COMBINE_DRF_OUTPUT.out.low_mapq_out, ref, ref_index,
                     PREPARE_ANNOTATION_BED.out.prepped_anno_bed, sequencer,
                     ref_tag)
}
