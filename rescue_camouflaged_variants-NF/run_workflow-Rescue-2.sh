#!/bin/bash

output_name="Ensembl_hg38_release_105_ploidy_8"

nextflow run rescue_camo_variants.nf \
	-with-report "${output_name}"-report.html \
	-with-trace "${output_name}"-trace.txt \
	-with-timeline "${output_name}"-timeline.html \
	-resume \
	--current_ref_fasta "/project/mteb223_uksr/sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
	--extraction_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla_illuminaRL100_Original_ADSP_samples-2022_06_06-17.09.25/05-CREATE_BED_FILE/illuminaRL100.1KGenomes_hg38_2015-GRCh38_full_analysis_set_plus_decoy_hla.camo.extraction.sorted.bed" \
	--masked_ref_fasta "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_Original_10_ADSP_samples-2022_06_05/06-MASK_GENOME/illuminaRL100.Ensembl_hg38_release_105.fa" \
	--unmasked_ref_fasta "/project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_105/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
	--masked_ref_tag "Ensembl_hg38_release_105_Original_10_ADSP_samples" \
	--input_sample_path "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/rescue_camouflaged_variants-NF/test_data/SMN_variants_samples" \
	--sample_input_tag "SMN_variants_samples_ALL" \
	--mask_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_Original_10_ADSP_samples-2022_06_05/05-CREATE_BED_FILE/illuminaRL100.Ensembl_hg38_release_105.camo.mask_bed.sorted.bed" \
	--gatk_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_Original_10_ADSP_samples-2022_06_05/05-CREATE_BED_FILE/illuminaRL100.Ensembl_hg38_release_105.camo.GATK.all_camo_regions.bed" \
	--annotation_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_Original_10_ADSP_samples-2022_06_05/04-PREPARE_ANNOTATION_BED/Homo_sapiens.GRCh38.105.annotation.bed" \
	--camo_annotations "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/Ensembl_hg38_release_105_Original_10_ADSP_samples-2022_06_05/05-CREATE_BED_FILE/illuminaRL100.Ensembl_hg38_release_105.camo_annotations.txt" \
	--n_samples_per_batch 50 \
	--rescue_gene_elements "ALL" \
	--max_repeats_to_rescue 5 \
	--ploidy_to_use 8 \
	--clean_tmp_files "false" \
	--DRF_interval_length 5000000
	
