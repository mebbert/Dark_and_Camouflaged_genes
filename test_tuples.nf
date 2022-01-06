// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2

// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##

params.gvcf_result_dir = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/nextflow/test_data/gvcfs'
params.gatk_bed = ''
params.mask_ref_prefix = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/nextflow/test_data/results/hg38_camo_mask/hg38-camo_mask'
params.result_dir_main = '/pscratch/mteb223_uksr/rescue_camo_variants/maddy_testing_scripts/nextflow/test_data/results/SAMPLES_WES'

log.info """\
 TEST PIPELINE
 ==========================================
 """

// Import Modules
include {test_tuples} from './modules/test_tuples_file.nf'

// Define initial files and channels
gatk_bed = file(params.gatk_bed)


workflow{
	test_tuples(params.gvcf_result_dir, gatk_bed, params.result_dir_main, params.mask_ref_prefix)
	test_tuples.out.test_tuples_out.view()
}
