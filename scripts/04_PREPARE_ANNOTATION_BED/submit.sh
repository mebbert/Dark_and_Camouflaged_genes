#!/bin/bash

## Update Ensembl's GFF3 Gene Annotations
## to remove overlapping gene body elements 
## (I.E. collapse down multiple transcripts per gene down to one,
## remove exon labels from CDS/UTR labels, add Introns, and remove 
## antisense or overlapping genes)

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
