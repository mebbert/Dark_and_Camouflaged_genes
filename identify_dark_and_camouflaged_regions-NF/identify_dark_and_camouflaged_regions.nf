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
params.align_to_ref = "${projectDir}/../sequencing_resources/references/Ensembl/hg38_release_93/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
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
params.align_to_ref_tag = 'Ensembl_GRCh38_release_93'
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
// params.input_sample_path = "${projectDir}/original_ADSP_samples/*.cram"

/*
 * This string can be used to help name the results folder (if provided and included when defining
 * the results parameter. 
 */
params.sample_input_tag = "Test_samples"
// params.sample_input_tag = "Original_ADSP_samples"

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
params.mask_ref_prefix = 'hg38_camo_mask'

log.info """\
 IDENTIFY DARK AND CAMO REGIONS PIPELINE
 ==========================================
 original reference             : ${params.original_ref}
 align_to_ref                   : ${params.align_to_ref}
 align_to_ref_tag               : ${params.align_to_ref_tag}
 input_sample_path              : ${params.input_sample_path}
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


align_to_ref = file(params.align_to_ref)
original_ref = file(params.original_ref)
align_to_ref_tag = params.align_to_ref_tag
output_format = params.output_format
align_to_gff = params.align_to_gff
sequencer_tag = params.sequencer_tag
DRF_jar = file("${projectDir}/bin/DarkRegionFinder.jar")


/*
 * Collect sample sam/bam/cram files from file input channel
 */
def collect_sample_files( input_dir ) {

    /*
     * Only accept sam/bam/cram
     */ 
    def pattern = ~/.*(\.sam|\.bam|\.cram)/
    def results = []

    input_dir.eachFileMatch(pattern) { item ->
        results.add( item )
    }

    return results
}

def calculate_ref_genome_size( ref ) {
    def ref_faidx = file("${ref}.fai")
    // println "faidx: ${ref_faidx}"
    def line, contig_length
    def long total_length = 0
    def String[] toks
    ref_faidx.withReader { reader ->
        while ((line = reader.readLine()) != null) {
            // println "#########"
            // println "${line}"
            toks = line.split("\t")
            contig_length = toks[1].toInteger()
            total_length += contig_length

            // println "Contig length: ${contig_length}"
            // println "Cumulative length: ${total_length}"
            // println ""
        }
    }
    return total_length
}

def create_intervals( ref, interval_length ) {
    
    def ref_faidx = file("${ref}.fai")
    def tmp_interval, remaining_contig_length, intervals = []
    ref_faidx.withReader { reader ->
        while ((line = reader.readLine()) != null) {
            // println "#################"
            // println "${line}"
            toks = line.split("\t")
            contig = toks[0]
            contig_length = toks[1].toInteger()

            remaining_contig_length = contig_length

            /*
             * Split contig into intervals of size 'interval_length'
             */
            // println "Intervals:"
            while (remaining_contig_length > 0) {
                
                /*
                 * 0-based start & end
                 */
                start = 0 + (contig_length - remaining_contig_length)

                if( remaining_contig_length < interval_length ) {
                    end = contig_length
                }
                else {
                    end = start + interval_length
                }
                tmp_interval = "${contig}:${start}-${end}"
                // println tmp_interval
                intervals.add(tmp_interval)

                remaining_contig_length -= interval_length
            }
            // println ""
        }
    }

    return intervals

}


workflow{



    /*
     * Create channel from the path provided, allowing only sam/bam/cram files.
     */
    // input_sample_path = file( params.input_sample_path )
    // input_sample_path = file( "${projectDir}/test_data/ADSP_sample_crams/*.cram" )
    input_sample_path = file( "${projectDir}/test_data/ADSP_sample_crams/A-CUHS-CU003023-test.cram" )
    input_sample_file_ch = Channel.fromPath(input_sample_path, checkIfExists: true)
                        .filter( ~/.*(\.sam|\.bam|\.cram)/ )

    /*
     * Prefix for various output files.
     */
    file_prefix = "${sequencer_tag}.${align_to_ref_tag}"

    /*
     * Realign input sample files
     */
//	REALIGN_SAMPLES_PROC(input_sample_file_ch, align_to_ref, align_to_ref_tag, original_ref,
//                     output_format)
    reads_per_job = 100_000
	REALIGN_SAMPLES_WF(input_sample_file_ch, reads_per_job, original_ref, align_to_ref, align_to_ref_tag,
                       output_format)

    /*
     * Run 'Dark Region Finder' to create summary statistics for every position in the
     * 'align_to_ref'. This process is run  split into a different process for each input sample
     * file.
     *
     * DRF is begging to be parallelized, but we haven't done that. Very slow as it is.
     */

    /* Determine total length of genome */
//    align_to_ref_length = calculate_ref_genome_size( align_to_ref )
//    n_DRF_jobs = 400
//    interval_length = align_to_ref_length.intdiv(n_DRF_jobs)

    // println "interval length: ${interval_length}"
    // intervals = Channel.from( create_intervals( align_to_ref, interval_length ) )

    interval_length = 150000000

 	// RUN_DRF(REALIGN_SAMPLES.out.final_alignment, align_to_ref, DRF_jar, intervals)
 	// RUN_DRF("${projectDir}/./results/Ensembl_GRCh38_release_93_illuminaRL100_Original_ADSP_samples-2022_02_03-18.07.32/REALIGN_SAMPLES/A-CUHS-CU000208-BL-COL-56227BL1.Ensembl_GRCh38_release_93.bam", align_to_ref, DRF_jar, intervals)
 	// RUN_DRF_WF("${projectDir}/./results/Ensembl_GRCh38_release_93_illuminaRL100_Original_ADSP_samples-2022_02_03-18.07.32/REALIGN_SAMPLES/A-CUHS-CU000208-BL-COL-56227BL1.Ensembl_GRCh38_release_93.bam", align_to_ref, DRF_jar, interval_length)

    // samples_ch = Channel.fromPath("${projectDir}/./results/Ensembl_GRCh38_release_93_illuminaRL100_Original_ADSP_samples-2022_02_03-18.07.32/REALIGN_SAMPLES/*.bam")
 	// RUN_DRF_WF(samples_ch, align_to_ref, DRF_jar, interval_length)

// 
//     /*
//      * Combine DRF output from all samples while also separating the DRF results into separate .bed
//      * files for 'dark-by-depth' and 'dark-by-MAPQ'. This w across all samples run previously. The
//      * string passed in is used for naming final results output.
//      */
// 	COMBINE_DRF_OUTPUT_PROC(RUN_DRF.out.low_mapq_bed.collect(), file_prefix)
// 
//     /*
//      * Calculate bam/cram metrics
//      */
// //	step_03_CALC_BAM_METRIC()
// 
//     /*
//      * Prepare the input annotation .gff3 file for use with defining camouflaged regions.
//      */
// 	PREPARE_ANNOTATION_BED_PROC(align_to_gff)
// 
//     /*
//      * Create final output files defining 'dark' and 'camouflaged' regions, along with various
//      * statistics.
//      */
// 	CREATE_BED_FILE_PROC(COMBINE_DRF_OUTPUT.out.low_depth_out,
//                      COMBINE_DRF_OUTPUT.out.low_mapq_out, align_to_ref, 
//                      PREPARE_ANNOTATION_BED.out.prepped_anno_bed, sequencer_tag,
//                      align_to_ref_tag)
// 
//     /*
//      * Mask the 'align_to_ref' for use in rescuing camouflaged variants from camouflaged regions in
//      * existing short-read sequencing data sets.
//      */
//  MASK_GENOME_PROC(CREATE_BED_FILE.out.align_to_bed, align_to_ref, file_prefix)

}
