#!/bin/bash

## 00_GET_BAMS
##
## After downloading the bams / fastqs files described in our
## Manuscript this submit will script will realign the data to
## get it in the correct format and aligned to the correct builds
## to use down in the camo analysis

## bwa mem is used to realign all the illumina short reads to 
## the three different references: b37, hg38, and hg38+alt

## minimap2 is used to realign PacBio and ONT long reads to 
## the correct references. ONT is also downsampled to ~50x coverage

## longranger is used to align 10X fastqs to b37 and hg38 reference


####################################################################
## TODO: Fill in correct paths before running                     ##
####################################################################
il100_B37_Dir="ADSP_WGS_BAMS/" #directory containing the 10 illumina100 sample bams downloaded from ADSP (aligned to b37/hg19 by default)
il250_Hg38_Alt_Dir="1kGP_HIGH_COVERAGE_BAMS/" #Directory contianing the 10 1000Genome 250 readlength bams (aligned by default to hg38+alt)

TenX_Fastq_Dirs="HG00512_10X_Fastqs" #Dir containing all the HG00512 Fastqs downloaded from 10X Genomics support page

PacBio_Bam="HG005_PacBio.hg38.bam" #HG005 bam downloaded from NIST GIAB (aligned to hg38 by default)

ONT_Fastq_Dir="Cliveome3_fastqs/" #Fastqs from Cliveome3 promethion downloaded from ONT github

b37_ref="Homo_sapiens.GRCh37.fa" # b37 reference
hg38_ref="Homo_sapiens.GRCh38.fa" # hg38 reference
hg38_alt_ref="GRCh38_full_analysis_set_plus_decoy_hla.fa" #hg38+alt reference 

TenX_b37_ref="refdata-b37-2.1.0/"
TenX_hg38_ref="refdata=GRCh38-2.1.0/"

#Where to output bams for each technology
IL100_OUT="../results/illuminaRL100/alignments"
IL250_OUT="../results/illuminaRL250/alignments"
TENX_OUT="../results/10X/alignments"
PACBIO_OUT="../results/PacBio/alignments"
ONT_OUT="../results/ONT/alignments"
####################################################################


#Realign Illumina RL 100 ADSP bams
for bam in $il100_B37_Dir/*.bam
do
	qsub realign_bwa.ogs $bam $hg38_ref "hg38" $IL100_OUT/hg38
	qsub realign_bwa.ogs $bam $hg38_alt_ref "hg38+alt" $IL100_OUT/hg38_alt
done

#Realign Illumina RL 250 1KGP bams
for bam in $il250_Hg38_Alt_Dir/*.bam
do
	qsub realign_bwa.ogs $bam $b37_ref "b37" $IL250_OUT/b37
	qsub realign_bwa.ogs $bam $hg38_ref "hg38" $IL250_OUT/hg38
done

#Align 10X Fastq to refs using longranger
qsub run_longranger.ogs $TenX_Fastq_Dirs $TenX_b37_ref $TENX_OUT/b37
qsub run_longranger.ogs $TenX_Fastq_Dirs $TenX_hg38_ref $TENX_OUT/hg38

#Realign PacBio Bam to other genome builds
realign_minimap2.ogs $PacBio_Bam $b37_ref "b37" $PACBIO_OUT/b37
realign_minimap2.ogs $PacBio_Bam $hg38_alt_ref "hg38+alt" $PACBIO_OUT/hg38+alt

# Align Cliveome3 fastqs to all three builds
align_ONT_minimap2.ogs $ONT_Fastq_Dir $b37_ref "b37" $ONT_OUT/b37
align_ONT_minimap2.ogs $ONT_Fastq_Dir $hg38_ref "hg38" $ONT_OUT/hg38
align_ONT_minimap2.ogs $ONT_Fastq_Dir $hg38_alt_ref "hg38+alt" $ONT_OUT/hg38+alt
