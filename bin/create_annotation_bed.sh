#!/bin/bash

## Input a GFF3 gene annotation file downloaded from Ensembl 
# computes genome arithmetic on gff3 to return a bed where overlapping gene body elements are
# removed (transcripts are collapsed down and overlapping genes and antisense genes removed)

GFF3=$1
TMP="tmp/annotation/"

mkdir -p $TMP

base="${GFF3##*/}"
extension="${base##*.}"
if [[ $extension != "gff3" ]]
then
	echo "ERROR: Input File must be a gff3 with a .gff3 extension in the file name"
	exit 1
fi

FILTERED_GFF3="${TMP}/${base//.gff3/.filtered.gff3}"
echo "`date` awk '\$9 ~ \"^ID=gene:\"' $GFF3 | bedtools cluster -s | filter_gff3.py $GFF3 > $FILTERED_GFF3"
awk '$9 ~ "^ID=gene:"' $GFF3 | \
	bedtools cluster -s | \
	filter_gff3.py $GFF3 |
	sort -k1,1 -k4,4n > $FILTERED_GFF3

genes="${TMP}/${base//.gff3/.justGenes.bed}"
echo "`date` awk '\$9 ~ \"^ID=gene:\"' $FILTERED_GFF3 > $genes"
awk '$9 ~ "^ID=gene:" {print $1"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7"\t"$9}' $FILTERED_GFF3 > $genes

CDS="${TMP}/${base//.gff3/.cds.gff3}"
echo "`date` awk '\$3 == \"CDS\"' $FILTERED_GFF3 | sort -k1,1 -k4,4n > $CDS"
awk '$3 == "CDS"' $FILTERED_GFF3 | \
	sort -k1,1 -k4,4n > $CDS

UTRs="${TMP}/${base//.gff3/.utrs.bed}"
echo "`date` awk '\$3 ~ \"UTR\"' $FILTERED_GFF3 | sort -k1,1 -k4,4n |  bedtools merge -c 3,6,7 -o distinct  > $UTRs"
awk '$3 ~ "UTR"' $FILTERED_GFF3 | \
	sort -k1,1 -k4,4n | \
	bedtools merge -c 3,6,7 -o distinct | \
	bedtools subtract \
		-a - \
		-b $CDS \
		> $UTRs

exons="${TMP}/${base//.gff3/.justExons.bed}"
echo "`date` awk '$\3 == \"exon\"' $FILTERED_GFF3 |  sort -k1,1 -k4,4n |  bedtools merge -c 3,6,7 -o distinct |  bedtools subtract  -a -  -b $UTRs  |  cat - $UTRs |  sort -k1,1 -k2,2n -k3,3n > $exons"
awk '$3 == "exon"' $FILTERED_GFF3 | \
	sort -k1,1 -k4,4n | \
	bedtools merge -c 3,6,7 -o distinct | \
	bedtools subtract \
		-a - \
		-b $UTRs | \
	cat - $UTRs | \
	sort -k1,1 -k2,2n -k3,3n > $exons

annotation_bed=${base//.gff3/.annotation.bed}
echo "`date` bedtools intersect  -a $genes -b $exons -loj > test.txt python prepare_gene_level_gtf.py  > $annotation_bed"

bedtools intersect  -a $genes -b $exons -s -loj | \
	sort -k1,1 -k2,2n -k9,9n | \
	prepare_annotation_bed.py > $annotation_bed
