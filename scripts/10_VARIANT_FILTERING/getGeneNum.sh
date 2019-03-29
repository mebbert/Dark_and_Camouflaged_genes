FILTERED_VCF=$1
ANNOTATION=$2
OUT_DIR=$3

awk '$6 == "PASS" && $9 >= 2.0' $FILTERED_VCF | \
	bedtools intersect \
		-a $ANNOTATION \
		-b - | \
	awk '$4 == "CDS" {print $NF}' | \
	python getGeneNumber.py \
	> ${OUT_DIR}/geneNumber.QD2.0.txt &

awk '$6 == "PASS" && $9 >= 3.0' $FILTERED_VCF | \
	bedtools intersect \
		-a $ANNOTATION \
		-b - | \
	awk '$4 == "CDS" {print $NF}' | \
	python getGeneNumber.py \
	> ${OUT_DIR}/geneNumber.QD3.0.txt
