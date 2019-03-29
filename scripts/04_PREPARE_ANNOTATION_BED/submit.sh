#!/bin/bash

## 04_PREPARE_ANNOTATION_BED
## Step 4 of our analysis converts Ensembl GFF3 Gene Annotations
## from a transcript level to a gene level, using genome arithmetic
## We only consider normal transcripts, ignoring all nonsense-mediated decay
## or retained intron transcripts when calculating gene-level annotations

## We also format annotation beds, removing overlapping gene body elements,
## remove exon labels from CDS/UTR labels, add Intron annotations, and remove 
## antisense or overlapping genes (i.e. miRNA genes present within the
## intron of a larger gene)

## We run these scripts for both b37 and hg38. For the hg38 annotations we
## must then add "chr" to the front of each chrom so it is consistent with
## hg38 contig naming.


###################################
## TODO: Fill in Correct Paths   ##
###################################
B37_GFF3="Homo_sapiens.GRCh37.87.gff3" #Path to Ensembl HG37 gene annotation GFF3
B37_RESULT_DIR="../../results/annotation/b37" #Dir to write resulting bed

B38_GFF3="Homo_sapiens.GRCh38.93.gff3" #Path to Ensembl HG38 gene annotation GFF3
B38_RESULT_DIR="../../results/annotation/b38"
####################################

bash create_annotation_bed.sh $B37_GFF3 $B37_RESULT_DIR
bash create_annotation_bed.sh $B38_GFF3 $B38_RESULT_DIR
