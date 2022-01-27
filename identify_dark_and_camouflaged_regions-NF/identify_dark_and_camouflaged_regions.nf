// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2

// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --ref reference.fa) ##

/*
 * Path to the reference the sample(s) are *currently* aligned to. This is only used if the input
 * sample files are .cram because bedtools bamtofastq needs the reference to convert the .cram to
 * fastq for re-alignment.
 *
 * If the input sample files are in .sam or .bam format, this paramater can be set to empty without
 * breaking. Leaving original_ref empty when input files are in .cram format is untested and could
 * lead to erroneous results.
 *
 * IMPORTANT: All input samples must already be aligned to the same reference.
 */
params.original_ref = "${projectDir}/../sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" //add a default value

/*
 * This is the reference genome that samples will be re-aligned to and that will
 * be used to define 'dark' and 'camouflaged' regions.
 *
 * IMPORTANT: It is essential to be very specific about the reference genome that is being used to
 * define 'dark' and 'camouflaged' regions. There are many different versions--even within the
 * "same" build (e.g., GRCh38)--that will have wildly different results. Even seemingly
 * insignificant modifications can have a major impact on which regions will be dark or camouflaged.
 * Specifically, including alternate contigs will have a major impact. See our
 * paper for details.
 */
params.align_to_ref = "${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

/*
 * This string is used in various output files to help identify
 * exactly which genome version is being used to define 'dark' and 'camouflaged' regions.
 *
 * IMPORTANT: We strongly recommend this tag be very specific. i.e., not just 'hg38', but which
 * specific version of 'hg38', including where it was prepared (e.g., 'Ensembl_GRCh38_release_93').
 * Even seemingly insignificant modifications can have a major impact on which regions will be
 * dark or camouflaged. Specifically, including alternate contigs will have a major impact. See our
 * paper for details.
 */
params.align_to_ref_tag = 'Ensembl_GRCh38_release_93'

/*
 * Path the the folder where input sample files are located. These files must be in either .sam,
 * .bam, or .cram format (or any combination thereof). 
 * 
 * The input can simply a path to the folder (e.g., /path/to/folder/; w/ or w/o ending '/'), or it
 * can specify a glob to limit the intput to specific files (e.g., /path/to/folder/*.cram,
 * /path/to/folder/<specific_file>.bam, etc.).
 */
// params.input_sample_folder = "${projectDir}/../samples/ADSP/input_sample_folder/*.cram"
// params.input_sample_folder = "${projectDir}/test_data/ADSP_sample_input_sample_folder/*"
params.input_sample_folder = "${projectDir}/original_ADSP_samples/"

/*
 * Specifies the final output format for re-aligned samples. Can be either '.bam' or '.cram' (w/ or
 * w/o the '.'
 */
params.output_format = 'bam'

/*
 * Path to the .gff3 file with all gene and transcript annotations. This MUST be specfic to the
 * 'align_to_ref' that input samples will be re-aligned to, otherwise the results will be entirely
 * erroneous. 
 */
params.align_to_gff = "${projectDir}/../sequencing_resources/annotations/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.93.gff3"

/*
 * Specifies the sequencer/technology that the input files were generated on (e.g.,
 * 'Illumina_RL_100' to specify that these samples were sequenced on Illumina technology with 100
 * nucleotide read lengths (i.e., 100 cycles).
 * 
 * Similar to the 'align_to_ref_tag', this is used in naming various output files to clearly
 * identify how the 'dark' and 'camouflaged' regions are defined. Read length is as important as
 * which reference genome is used to define 'dark' and 'camouflaged' regions!
 */
params.sequencer_tag = 'illuminaRL100'

/*
 * Path where final results should be saved. The workflow will generate subfolders for each of the
 * steps below with their respective outputs.
 *
 * IMPORTANT: NextFlow is picky about relative paths. i.e., simply using 'results' is not allowed.
 * It must be either an absolute path or a relative path that includes a preceing './' or the global
 * variable "${project_dir}" (e.g., '${project_dir}/').
 */
params.results_dir = './results'

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
 align_to_ref_tag                : ${params.align_to_ref_tag}
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
align_to_ref_tag = params.align_to_ref_tag
sequencer_tag = params.sequencer_tag
run_drf_jar = file("${projectDir}/bin/DarkRegionFinder.jar")


workflow{

    /*
     * Prefix for various output files.
     */
    file_prefix = "${sequencer_tag}.${align_to_ref_tag}"

    /*
     * Realign input sample files
     */
	REALIGN_SAMPLES(samples, align_to_ref, align_to_ref_tag,
                    original_ref, output_format)

    /*
     * Run 'Dark Region Finder' to create summary statistics for every position in the
     * 'align_to_ref'. This process is run  split into a different process for each input sample
     * file.
     *
     * DRF is begging to be parallelized, but we haven't done that. Very slow as it is.
     */
	RUN_DRF(REALIGN_SAMPLES.out.final_alignments, align_to_ref, run_drf_jar)

    /*
     * Combine DRF output across all samples run previously. The string passed in is used for naming
     * final results output.
     */
	COMBINE_DRF_OUTPUT(RUN_DRF.out.low_mapq_bed.collect(), file_prefix)

    /*
     * Calculate bam/cram metrics
     */
//	step_03_CALC_BAM_METRIC()

    /*
     * Prepare the input annotation .gff3 file for use with defining camouflaged regions.
     */
	PREPARE_ANNOTATION_BED(params.align_to_gff)

    /*
     * Create final output files defining 'dark' and 'camouflaged' regions, along with various
     * statistics.
     */
	CREATE_BED_FILE(COMBINE_DRF_OUTPUT.out.low_depth_out,
                     COMBINE_DRF_OUTPUT.out.low_mapq_out, align_to_ref, ref_index,
                     PREPARE_ANNOTATION_BED.out.prepped_anno_bed, sequencer_tag,
                     align_to_ref_tag)

    /*
     * Maske the 'align_to_ref' for use in rescuing camouflaged variants from camouflaged regions in
     * existing short-read sequencing data sets.
     */
    MASK_GENOME(CREATE_BED_FILE.out.align_to_bed, align_to_ref, file_prefix)

}
