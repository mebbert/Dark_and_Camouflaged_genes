#!/bin/bash

output_name="T2T_CHM13_v2.0"

nextflow run rescue_camo_variants.nf \
	-with-report "${output_name}"-report.html \
	-with-trace "${output_name}"-trace.txt \
	-with-timeline "${output_name}"-timeline.html \
	-resume \
	--current_ref_fasta "/project/mteb223_uksr/sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
	--extraction_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/T2T_CHM13_v2.0_illuminaRL100_Original_ADSP_samples-2022_06_01-16.34.59/05-CREATE_BED_FILE/illuminaRL100.T2T_CHM13_v2.0.camo.extraction.sorted-SMN_only.bed" \
	--masked_ref_fasta "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/T2T_CHM13_v2.0_illuminaRL100_Original_ADSP_samples-2022_06_01-16.34.59//06-MASK_GENOME/illuminaRL100.T2T_CHM13_v2.0.fa" \
	--unmasked_ref_fasta "/project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_105/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
	--masked_ref_tag "T2T_CHM13_v2.0.Original_10_ADSP_samples" \
	--input_sample_path "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/rescue_camouflaged_variants-NF/test_data/SMN_variants_samples" \
	--sample_input_tag "SMN_samples_ADSP" \
	--mask_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/T2T_CHM13_v2.0_illuminaRL100_Original_ADSP_samples-2022_06_01-16.34.59/05-CREATE_BED_FILE/illuminaRL100.T2T_CHM13_v2.0.camo.mask_bed.sorted-SMN_only.bed" \
	--gatk_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/T2T_CHM13_v2.0_illuminaRL100_Original_ADSP_samples-2022_06_01-16.34.59/05-CREATE_BED_FILE/illuminaRL100.T2T_CHM13_v2.0.camo.GATK.CDS_regions_only-SMN_only.bed" \
	--annotation_bed "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/T2T_CHM13_v2.0_illuminaRL100_Original_ADSP_samples-2022_06_01-16.34.59/04-PREPARE_ANNOTATION_BED/CHM13.v2.0.sorted_attributes.annotation-SMN_only.bed" \
	--camo_annotations "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/results/T2T_CHM13_v2.0_illuminaRL100_Original_ADSP_samples-2022_06_01-16.34.59/05-CREATE_BED_FILE/illuminaRL100.T2T_CHM13_v2.0.camo_annotations-SMN_only.txt" \
	--n_samples_per_batch 50 \
	--rescue_gene_elements "CDS" \
	--max_repeats_to_rescue 5 \
	--clean_tmp_files "false" \
	--DRF_interval_length 5000000
	
