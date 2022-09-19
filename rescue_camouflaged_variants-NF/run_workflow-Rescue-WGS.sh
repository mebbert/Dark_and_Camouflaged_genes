#!/bin/bash
#SBATCH --time=60:15:00                                         # Time limit for the job (REQUIRED).
#SBATCH --job-name=NF_Parent_Rescue_Camo_Vars      # Job name
#SBATCH --ntasks=1                                              # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                                                # Total memory requested
#SBATCH --partition=normal                                      # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err                                         # Error file for this job.
#SBATCH -o slurm-%j.out                                         # Output file for this job.
#SBATCH -A coa_mteb223_uksr                                     # Project allocation account name (REQUIRED)

module load ccs/java/jdk1.8.0_202

output_name="Ensembl_hg38_release_105_ploidy_8_WGS"


	#-with-report "${output_name}"-report.html \
	#-with-trace "${output_name}"-trace.txt \
	#-with-timeline "${output_name}"-timeline.html \

nextflow run rescue_camo_variants.nf \
	-resume \
	--current_ref_fasta "/project/mteb223_uksr/sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
	--extraction_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_09_06-09.52.01/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo.extraction.sorted.bed" \
	--masked_ref_fasta "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_illuminaRL100_Original_10_ADSP_samples-2022_09_05-09.51.20/06-MASK_GENOME/illuminaRL100.Ensembl_hg38_release_105.fa" \
	--unmasked_ref_fasta "/project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_105/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
	--masked_ref_tag "Ensembl_hg38_release_105_Original_10_ADSP_samples" \
	--input_sample_path "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/rescue_camouflaged_variants-NF/test_data/WGS_crams" \
	--sample_input_tag "WGS_adsp_samples" \
	--mask_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_illuminaRL100_Original_10_ADSP_samples-2022_09_05-09.51.20/05-CREATE_BED_FILE/illuminaRL100.Ensembl_hg38_release_105.camo.mask_bed.sorted.bed" \
	--gatk_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_illuminaRL100_Original_10_ADSP_samples-2022_09_05-09.51.20/05-CREATE_BED_FILE/illuminaRL100.Ensembl_hg38_release_105.camo.GATK.all_camo_regions.bed" \
	--annotation_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_illuminaRL100_Original_10_ADSP_samples-2022_09_05-09.51.20/04-PREPARE_ANNOTATION_BED/Homo_sapiens.GRCh38.105.annotation.bed" \
	--camo_annotations "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_illuminaRL100_Original_10_ADSP_samples-2022_09_05-09.51.20/05-CREATE_BED_FILE/illuminaRL100.Ensembl_hg38_release_105.camo_annotations.txt" \
	--n_samples_per_batch 50 \
	--rescue_gene_elements "ALL" \
	--max_repeats_to_rescue 5 \
	--ploidy_to_use 8 \
	--clean_tmp_files "false" \
	--DRF_interval_length 5000000
	
