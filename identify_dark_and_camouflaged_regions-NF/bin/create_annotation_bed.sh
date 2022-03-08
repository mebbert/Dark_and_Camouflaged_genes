#!/bin/bash

# Enable bash debugging to log all commands
set -x

## Input a GFF3 gene annotation file downloaded from Ensembl 
# computes genome arithmetic on gff3 to return a bed where overlapping gene body elements are
# removed (transcripts are collapsed down and overlapping genes and antisense genes removed)

GFF3=$1
TMP="annotation/"

mkdir -p $TMP

base="${GFF3##*/}"
extension="${base##*.}"
if [[ $extension != "gff3" ]]
then
	echo "ERROR: Input File must be a gff3 with a .gff3 extension in the file name"
	exit 1
fi

# Create a filtered GFF file that contains only a single gene entry for any given gene.
# Essentially, because the GFF file contains gene definitions from multiple consortia
# (e.g., havana, ensembl), there are multiple and differing gene definitions in the file.
# We chose to select the largest of the gene definitions for any given gene (based on start
# and end. This will also include all transcripts from that gene definition.
FILTERED_GFF3="${TMP}/${base//.gff3/.filtered.gff3}"
# echo "`date` awk '\$9 ~ \"^ID=gene:\"' $GFF3 | bedtools cluster -s | filter_gff3.py $GFF3 > $FILTERED_GFF3"
# TODO: Keep an eye on whether the added sort causes a problems for reference genomes
# where chromosomes are sorted numerically rather than ascii-betically. 
awk '$9 ~ "^ID=gene:"' $GFF3 | \
	bedtools cluster -s | \
	filter_gff3.py $GFF3 | \
	sort -k1,1 -k4,4n > $FILTERED_GFF3

# Isolate just the gene lines from the filtered GFF into a single
# genes GFF file.
genes="${TMP}/${base//.gff3/.justGenes.bed}"
# echo "`date` awk '\$9 ~ \"^ID=gene:\"' $FILTERED_GFF3 > $genes"
awk '$9 ~ "^ID=gene:" {print $1"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7"\t"$9}' $FILTERED_GFF3 > $genes

# Isolate just the CDS lines from the filtered GFF
CDS="${TMP}/${base//.gff3/.cds.gff3}"
# echo "`date` awk '\$3 == \"CDS\"' $FILTERED_GFF3 | sort -k1,1 -k4,4n > $CDS"
awk '$3 == "CDS"' $FILTERED_GFF3 | \
	sort -k1,1 -k4,4n > $CDS

# Isolate just the UTRs from the filtered GFF. In this case, however, we also
# subtract the CDS regions from the UTRs. i.e., because the CDS region is
# considered CDS in at least one transcript, we will always treat it as CDS
# and NOT a UTR.
UTRs="${TMP}/${base//.gff3/.utrs.bed}"
# echo "`date` awk '\$3 ~ \"UTR\"' $FILTERED_GFF3 | sort -k1,1 -k4,4n |  bedtools merge -c 3,6,7 -o distinct  > $UTRs"
awk '$3 ~ "UTR"' $FILTERED_GFF3 | \
	sort -k1,1 -k4,4n | \
	bedtools merge -s -c 3,6,7 -o distinct | \
	bedtools subtract \
		-a - \
		-b $CDS \
		-s > $UTRs

# Isolate just the exons (including UTRs) from the filtered GFF. In this case,
# we give preference to UTRs since CDS was given preference in the previous
# command. i.e., If we want to be able to look at all UTRs, this is all
# inclusive.
exons="${TMP}/${base//.gff3/.justExons.bed}"
# echo "`date` awk '$\3 == \"exon\"' $FILTERED_GFF3 |  sort -k1,1 -k4,4n |  bedtools merge -c 3,6,7 -o distinct |  bedtools subtract  -a -  -b $UTRs  |  cat - $UTRs |  sort -k1,1 -k2,2n -k3,3n > $exons"
awk '$3 == "exon"' $FILTERED_GFF3 | \
	sort -k1,1 -k4,4n | \
	bedtools merge -s -c 3,6,7 -o distinct | \
	bedtools subtract -s \
		-a - \
		-b $UTRs | \
	cat - $UTRs | \
	sort -k1,1 -k2,2n -k3,3n > $exons

annotation_bed=${base//.gff3/.annotation.bed}
# echo "`date` bedtools intersect  -a $genes -b $exons -loj > test.txt python prepare_gene_level_gtf.py  > $annotation_bed"

# Create a single GFF file that contains individual entries for every
# exon for a given gene, including all gene annotations. These annotations
# are not included in the exons file. 
bedtools intersect  -a $genes -b $exons -s -loj | \
	sort -k1,1 -k2,2n -k9,9n | \
	prepare_annotation_bed.py > $annotation_bed

# Check that $annotation_bed is not empty. Fail, if so.
anno_lines=$(cat $annotation_bed | wc -l)
if [[ "${anno_lines}" == 0 ]]; then
	echo "ERROR: $annotation_bed is empty!"
	exit 1
fi
