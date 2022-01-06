// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2

// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##

params.ref='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.ref_index='/project/mteb223_uksr/sequencing_resources/references/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
//params.align_to_bed = '/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/work/8e/5408b4aa7f4dd3264d56d0fd85afac/results_dir/illuminaRL100.hg38.camo.align_to.sorted.bed'
params.align_to_bed = '/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/test_data/small_align_to.bed'
params.realign_bed = '/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/test_data/illuminaRL100.b37.camo.realign.sorted.bed'
//params.realign_bed = '/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/work/8e/5408b4aa7f4dd3264d56d0fd85afac/results_dir/illuminaRL100.hg38.camo.realign.sorted.bed'
params.gatk_bed = '/mnt/gpfs3_amd/condo/mteb223/rescue_camo_variants/nextflow/test_data/illuminaRL100.b37.camo.GATK.bed'
params.ref_tag = 'hg38'
params.mask_ref_prefix = 'hg38_camo_mask'
params.filtered_bam_list = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/nextflow/test_data/ADSP_sample_crams/*'

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
include {step_06_MASK_GENOME} from './modules/06_MASK_GENOME.nf'
include {step_07_PREP; step_07_RUN_SAMPLES} from './modules/07_RUN_SAMPLES.nf'

// Define initial files and channels
ref = file(params.ref)
ref_idx = file(params.ref_index)
align_to_bed = file(params.align_to_bed)
realign_bed = file(params.realign_bed)
filtered_bams_list = params.filtered_bam_list
gatk_bed = params.gatk_bed


workflow{
	step_06_MASK_GENOME(align_to_bed, ref, ref_idx, params.mask_ref_prefix)
	step_07_PREP(filtered_bams_list)
	step_07_RUN_SAMPLES(params.realign_bed, params.gatk_bed, params.mask_ref_prefix, step_07_PREP.out.filtered_lists.flatten())

	
}
