#!/bin/bash

output_name="NCBI_GRCh38_no_alt_plus_hs38d1"

nextflow run identify_dark_and_camouflaged_regions.nf \
	-with-report "${output_name}"-report.html \
	-with-trace "${output_name}"-trace.txt \
	-with-timeline "${output_name}"-timeline.html \
	-resume \
	--original_ref "/project/mteb223_uksr/sequencing_resources/references/1KGenomes_hg38-2015/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
	--align_to_ref "/project/mteb223_uksr/sequencing_resources/references/NCBI/GRCh38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna" \
	--align_to_ref_tag "NCBI_GRCh38_no_alt_plus_hs38d1" \
	--input_sample_path "/mnt/gpfs3_amd/condo/mteb223/mlpa241/Dark_and_Camouflaged_genes/identify_dark_and_camouflaged_regions-NF/original_ADSP_samples/*.cram" \
	--sample_input_tag "Original_10_ADSP_samples" \
	--reads_per_run 50000000 \
	--output_format "bam" \
	--align_to_gff "/project/mteb223_uksr/sequencing_resources/annotations/Ensembl/hg38_release_105/Homo_sapiens.GRCh38.105.chr_added.gff3" \
	--sequencer_tag "illuminaRL100" \
	--mask_ref_prefix "NCBI_GRCh38_no_alt_plus_hs38d1" \
	--DRF_interva_length 5000000 
