import java.text.SimpleDateFormat
import java.io.File

/*
 * Make this pipeline a nextflow 2 implementation
 */
nextflow.enable.dsl=2

/*
 * Get current date and time
 */
def date = new Date()
def sdf = new SimpleDateFormat("yyyy_MM_dd-HH.mm.ss")
def time_stamp = sdf.format(date)

/*
 * Pipeline parameter default values, can be modified by user when calling pipeline on command line
 * (e.g. --ref reference.fa) ##
 */

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
// params.align_to_ref = "${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.align_to_ref = "/project/mteb223_uksr/sequencing_resources/references/Tanner_Original_Camo_References/refdata-GRCh38-2.1.0/fasta/genome.fa" \
// params.align_to_ref = "${projectDir}/../sequencing_resources/references/Gencode/release_31/GRCh38.primary_assembly/GRCh38.primary_assembly.genome.fa"

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
// params.align_to_ref_tag = 'Ensembl_GRCh38_release_93'
params.align_to_ref_tag = 'Tanner_Original_Camo_Reference'
// params.align_to_ref_tag = 'Gencode_release_31-primary_assembly'

/*
 * Path to input sample files including a 'wildcard' for which files to include. It cannot only be
 * the path to the folder. These files must be
 * in either .sam, .bam, or .cram format (or any combination thereof).
 * 
 * The glob must specify the intput (e.g., /path/to/folder/*, /path/to/folder/*.cram,
 * /path/to/folder/<specific_file>.bam, etc.). The input will be limited to .sam, .bam, and .cram
 * files before proceeding.
 */
params.input_sample_path = "${projectDir}/test_data/ADSP_sample_crams/*"
// params.input_sample_path = file( "${projectDir}/test_data/ADSP_sample_crams/A-CUHS-CU003023-test.cram" )
// params.input_sample_path = "${projectDir}/original_ADSP_samples/*.cram"

/*
 * This string can be used to help name the results folder (if provided and included when defining
 * the results parameter. 
 */
params.sample_input_tag = "Test_samples"
// params.sample_input_tag = "Original_ADSP_samples"

/*
 * Defines the number of reads per BWA alignment job. This number must be
 * divisible by 4 because reads in a .fastq file use 4 lines.
 *
 * NOTE: If this number results in a single .fastq file for the sample, Nextflow
 * will error out with something like the following:
 * Invalid method invocation `groupKey` with arguments: A-CUHS-CU003023-test (java.lang.String), 26766997 (java.lang.Long) on Nextflow type
 *
 * This error happens because `fastq_files.size()` returns the file size for the single file
 * returned, rather than reporting the size of the list of files, which is what is expected.
 * Nextflow channels seem inherently buggy to me.
 */
params.reads_per_run = 100_000
// params.reads_per_run = 50_000_000

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
 * variable "${projectDir}" (e.g., '${projectDir}/').
 */
params.results_dir = "./results/${params.align_to_ref_tag}_${params.sequencer_tag}_${params.sample_input_tag}-${time_stamp}"

/*
 * Prefix for masked reference files.
 */
params.mask_ref_prefix = params.align_to_ref_tag

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
 IDENTIFY DARK AND CAMO REGIONS PIPELINE
 ==========================================
 original reference             : ${params.original_ref}
 align_to_ref                   : ${params.align_to_ref}
 align_to_ref_tag               : ${params.align_to_ref_tag}
 input_sample_path              : ${params.input_sample_path}
 reads_per_run                  : ${params.reads_per_run}
 output_format                  : ${params.output_format}
 align_to_gff                   : ${params.align_to_gff}
 sequencer_tag                  : ${params.sequencer_tag}
 results_dir                    : ${params.results_dir}
 mask_ref_prefix                : ${params.mask_ref_prefix}
 """


/*
 * Import Modules
 */
include {REALIGN_SAMPLES_WF} from './modules/01-REALIGN_SAMPLES.nf'
include {RUN_DRF_WF} from './modules/02-RUN_DRF.nf'
include {COMBINE_DRF_OUTPUT_PROC} from './modules/03-COMBINE_DRF_OUTPUT.nf'
//include {step_03_CALC_BAM_METRICS} from './modules/03_CALC_BAM_METRICS.nf'
include {PREPARE_ANNOTATION_BED_PROC} from './modules/04-PREPARE_ANNOTATION_BED.nf'
include {CREATE_BED_FILE_PROC} from './modules/05-CREATE_BED_FILE.nf'
include {MASK_GENOME_PROC} from './modules/06-MASK_GENOME.nf'



/*
 * Let the magic begin.
 */
workflow{

    /* test that params.reads_per_run is divisible by 4 */
    if(params.reads_per_run  % 4 != 0){
        throw Exception("ERROR: reads_per_run must be divisible by 4. Please" +
            " specify a number that is divisible by 4.")
    }

    /*
     * Prefix for various output files.
     */
    file_prefix = "${params.sequencer_tag}.${params.align_to_ref_tag}"

    /*
     * Step 01: Realign input sample files
     */
	REALIGN_SAMPLES_WF(params.input_sample_path)

    /*
     * Step 02: Run 'Dark Region Finder' to create summary statistics for every
     * position in the 'align_to_ref'. This process is split into a different
     * process for each input sample file and genome intervals of size
     * 'params.DRF_interval_length'.
     */
 	RUN_DRF_WF(REALIGN_SAMPLES_WF.out, params.DRF_interval_length)
 

    /*
     * Step 03: Combine DRF output from all samples while also separating the
     * DRF results into separate .bed files for 'dark-by-depth' and
     * 'dark-by-MAPQ'. This will combine across all samples run previously.
     * The string passed in is used for naming final results output.
     */
    COMBINE_DRF_OUTPUT_PROC(RUN_DRF_WF.out.collect(), file_prefix)
 

    /*
     * Step 04: Prepare the input annotation .gff3 file for use with defining
     * camouflaged regions.
     */
 	PREPARE_ANNOTATION_BED_PROC(params.align_to_gff)
 
    /*
     * Step 05: Create final output files defining 'dark' and 'camouflaged'
     * regions, along with various statistics.
     */
 	CREATE_BED_FILE_PROC(COMBINE_DRF_OUTPUT_PROC.out.low_depth_out,
                      COMBINE_DRF_OUTPUT_PROC.out.low_mapq_out,
                      PREPARE_ANNOTATION_BED_PROC.out.prepped_anno_bed)
 
    /*
     * Step 06: Mask the 'align_to_ref' for use in rescuing camouflaged
     * variants from camouflaged regions in existing short-read sequencing data
     * sets.
     */
    MASK_GENOME_PROC(CREATE_BED_FILE_PROC.out.align_to_bed, file_prefix)

}
