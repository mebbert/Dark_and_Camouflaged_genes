#!/bin/bash

function usage {
	echo "usage: $ ./camo_gene_pipeline.sh [options...]"
	echo "Parameters:"
	echo "	-d, --low-depth-bed <FILE.bed>: DRF output bed file of low depth depth regions"
	echo " 	-m, --low-mapq-bed <FILE.bed>: DRF output bed file of low mapQ depth regions"
	echo "	-g, --genome_ref <FILE.fa>: path to reference genome fasta file"
	echo "	-a, --gene_annotations <FILE.bed>: path to a  gene annotation bed file (run create_gene_level_bed.sh on GFF3 file to create)"
	echo "	-s, --sequencer STR: sequencing technology used to generate bam files"
	echo "	-v, --genome_version STR: human genome version"
	echo "Options:"
	echo "	-t, --num_threads INT: number of threads to use (default = 1)"
	echo "	-r, --results DIR: directory where the results will be saved (default = Sequencer/GenomeVers)"
	exit 1
}

NUM_THREADS=1
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-d|--low-depth-bed)
			DEPTH_BED="$2"; shift; shift; ;;
		-m|--low-mapq-bed)
			MAPQ_BED="$2"; shift; shift; ;;
		-g|--genome)
			GENOME="$2"; shift; shift; ;;
		-a|--gene_annotations)
			ANNOTATION="$2"; shift; shift; ;;
		-s|--sequencer)
			SEQUENCER="$2"; shift; shift; ;;
		-v|--genome_version)
			GVERS="$2"; shift; shift; ;;
		-t|--num_threads)
			NUM_THREADS="$2"; shift; shift; ;;
		-r|--result_dir)
			RESULT_DIR="$2"; shift; shift; ;;
		*)
			echo "ERROR: unrecognized command line argument; exiting"
			usage
	esac
done

if [[ -z $DEPTH_BED || -z $MAPQ_BED || -z $GENOME || -z $ANNOTATION || -z $SEQUENCER || -z $GVERS ]]
then
	echo "ERROR: required arguments not specificied; exiting"
	usage
fi

if [[ ! -f "${GENOME}.fai" ]] 
then
	echo "ERROR: fasta index (.fai) must exist in same dir as genome reference"
	echo "run samtools faidx on $GENOME"
	usage
fi

if [[ -z $RESULT_DIR ]]
then
	RESULT_DIR=${GVERS}/${SEQUENCER}
fi

TMP_DIR="tmp/${SEQUENCER}.${GVERS}"
mkdir -p $TMP_DIR

##----Define Variables----------
faidx="${GENOME}.fai"
mapq_mass90="${TMP_DIR}/low_mapq-mass90.bed"
mapq_merged="${TMP_DIR}/${SEQUENCER}.${GVERS}.low_mapq-merged.bed"
depth_merged="${TMP_DIR}/${SEQUENCER}.${GVERS}.low_depth-merged.bed"
dark_merged="${TMP_DIR}/${SEQUENCER}.${GVERS}.dark-merged.bed"
mapq_lengths="${TMP_DIR}/${SEQUENCER}.${GVERS}.low_mapq.lengths.txt"
dark_loj="${TMP_DIR}/dark_loj.txt"
mapq_loj="${TMP_DIR}/mapq_loj.txt"
depth_loj="${TMP_DIR}/depth_loj.txt"
dark_annotations="${TMP_DIR}/${SEQUENCER}.${GVERS}.dark_annotations.txt"
mapq_annotations="${TMP_DIR}/${SEQUENCER}.${GVERS}.mapq_annotations.txt"
depth_annotations="${TMP_DIR}/${SEQUENCER}.${GVERS}.depth_annotations.txt"
camo_annotations="${TMP_DIR}/${SEQUENCER}.${GVERS}.camo_annotations.txt"
percent_dark="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_dark_genes.txt"
percent_mapq="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_low_mapq_genes.txt"
percent_depth="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_low_depth_genes.txt"
percent_camo="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_camo_genes.txt"
biotype_dark="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_dark_biotypes.txt"
biotype_mapq="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_mapq_biotypes.txt"
biotype_depth="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_depth_biotypes.txt"
biotype_camo="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_camo_biotypes.txt"
coding_dark="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_dark_coding_regions.txt"
coding_mapq="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_mapq_coding_regions.txt"
coding_depth="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_depth_coding_regions.txt"
coding_camo="${TMP_DIR}/${SEQUENCER}.${GVERS}.percent_camo_coding_regions.txt"
query="${TMP_DIR}/query/mapq_genes.query.fa"
blat_result="${TMP_DIR}/blat_result/blat.results.psl"
blat_log="${TMP_DIR}/blat_log/tmp.blat.log"
blat_bed=${blat_result//psl/bed}
mapped_blat_results="${TMP_DIR}/query.results.filtered.mapped.txt"
camo_bed="${TMP_DIR}/camo.bed"
align_to="${TMP_DIR}/align_to.bed"
realign="${TMP_DIR}/realign.bed"
camo_sorted="${TMP_DIR}/camo.sorted.bed"
alignto_sorted="${TMP_DIR}/${SEQUENCER}.${GVERS}.camo.align_to.sorted.bed"
realign_sorted="${TMP_DIR}/${SEQUENCER}.${GVERS}.camo.realign.sorted.bed"
gatk_bed="${TMP_DIR}/${SEQUENCER}.${GVERS}.camo.GATK.bed"
mapq_gene_bodies="${TMP_DIR}/low_mapq_gene_bodies.bed"
mapq_not_camo="${TMP_DIR}/${SEQUENCER}.${GVERS}.low_mapq.NOT_camo.bed"

#Merge coordinates for mapq and depth beds, removing regions that are less than 20bp long
echo "`date` bedtools merge -d 20-c 5 -o mean,median -i $DEPTH_BED > $depth_merged"
bedtools merge -d 20 -c 5 -o mean,median -i $DEPTH_BED | \
	remove_unassembled_contigs.py | \
	awk '{ if($3 - $2 > 20) print $0}' \
	> $depth_merged

echo "`date` bedtools merge -d 20 -c 5 -o mean,median -i $MAPQ_BED > $mapq_merged"
bedtools merge -d 20 -c 5 -o mean,median -i $MAPQ_BED | \
	remove_unassembled_contigs.py | \
	awk '{ if($3 - $2 > 20) print $0}' \
	> $mapq_merged

echo "`date` awk '{n +=1; print n,\$1,(\$3 - \$2)}' \$mapq_merged > \$mapq_lengths"
echo -e "region_id\tchrom\tlength" >> $mapq_lengths
awk '{n +=1; print n,$1,($3 - $2)}' $mapq_merged >> $mapq_lengths

cat $mapq_merged $depth_merged | \
	bedtools sort -i - -faidx $faidx | \
	bedtools merge -i - \
	> $dark_merged

## Intersect all dark regions with gene annotation be
echo "`date` bedtools intersect -a $ANNOTATION -b $dark_merged -loj | python annotate_regions.py > $dark_annotations"
bedtools intersect \
	-a $ANNOTATION \
	-b $dark_merged \
	-loj | \
	annotate_regions.py \
		$percent_dark $biotype_dark $coding_dark "dark" \
		> $dark_annotations

## INtersect the low depth dark regions with gene annotation bed
echo "`date` bedtools intersect -a $ANNOTATION -b $depth_merged -loj > $depth_loj | python annotate_regions.py > $depth_loj"
bedtools intersect \
	-a $ANNOTATION \
	-b $depth_merged \
	-loj | \
	annotate_regions.py \
		$percent_depth $biotype_depth $coding_depth "dark by depth" \
		> $depth_annotations

## Intersect the low mapq regions with gene annotation bed
echo "`date` bedtools intersect -a $ANNOTATION -b $mapq_merged -loj > $mapq_loj | python annotate_regions.py > $mapq_loj"
bedtools intersect \
	-a $ANNOTATION \
	-b $mapq_merged \
	-loj | \
	annotate_regions.py \
		$percent_mapq $biotype_mapq $coding_mapq "dark by MAPQ" \
		> $mapq_annotations

## Take mapq_file and create a blat query sequence list
mkdir -p "${TMP_DIR}/query"
cat $mapq_annotations | \
	grep -vE "^#" | \
	bedtools getfasta \
		-fi $GENOME \
		-bed - \
		-name \
		-fo $query

nLines=$(wc -l $query | awk '{print $1}')
stepSize=$(($nLines / $NUM_THREADS))
if [[ $(($stepSize % 2 )) == 1 ]]
then
	    stepSize=$(($stepSize + 1))
fi

rm -rf "${TMP_DIR}/blat_result"
mkdir -p "${TMP_DIR}/blat_result"
mkdir -p "${TMP_DIR}/blat_log"
for i in $(seq 1 $stepSize $nLines)
do
	sed "$((i)),$((i + $stepSize - 1))!d" $query > ${query}.${i}
	echo "`date` blat $GENOME ${query}.${i} -t=dna -q=dna blat.results.psl"
	blat $GENOME ${query}.${i} \
			-t=dna \
			-q=dna \
			-minIdentity=95 \
			-noHead \
			${blat_result}.${i} \
			&> ${blat_log}.${i} &
done
wait

## Combine the blat results together and format them
## into a bed file where score shows the sequence identity,
## removing blat hits that are less than 98% sequence identity
echo "`date` blatting complete: combining blat output"
cat ${blat_result}.* > $blat_result

echo "`date` awk -f score_blat_output.awk $blat_result > $blat_bed"
awk -f score_blat_output.awk \
	$blat_result \
	> $blat_bed

echo "`date` awk 'BEGIN{OFS="\t"}{ print $7,$8,$9,$4 }' $mapq_annotations > $mapq_gene_bodies"
awk -v "OFS=\t" '$1 !~ "^#" { print $7,$8,$9,$4 }' \
	$mapq_annotations \
	> $mapq_gene_bodies

## Intersect the blat output to original list of low_mapq regions
echo "`date` bedtools intersect -a $blat_bed -b $mapq_gene_bodies -loj | awk '$6 != "."' > $mapped_blat_results"
bedtools intersect \
	-a $blat_bed \
	-b $mapq_gene_bodies \
	-loj | \
	awk '$6 != "."' \
	> $mapped_blat_results

## From the scored intersected blat output, calculates the maps
## of which camo regions align to which other ones, and creates camo sets
# lists all the regions in  acamo set in the realign file and selects the one region form that set that we 
## will use to align to (written in the align_to file, all other regions in set will be masked)
echo "`date` python extract_camo_regions.py"
extract_camo_regions.py \
	$mapq_annotations \
	$mapped_blat_results \
	$realign \
	$align_to \
	$camo_bed

echo "`date` sorting camo bed files"

#extract_camo_regions prints out whole camo gene bodies,
#also same gene bodies may be printed twice, so first merge then
#intersect with merged_mapq to get just the Camo Region boundaries
bedtools sort -i $camo_bed -faidx $faidx | \
	bedtools merge -i - | \
	bedtools intersect -a - -b $mapq_merged > $camo_sorted
bedtools sort -i $align_to -faidx $faidx > $alignto_sorted
bedtools sort -i $realign  -faidx $faidx > $realign_sorted

## Create Camo Annotation table (intersecting camo regions to gene annotation bed)
echo "`date` bedtools intersect -a $ANNOTATION -b $camo_sorted -loj | python annotate_regions.py > $camo_annotations"
bedtools intersect \
	-a $ANNOTATION \
	-b $camo_sorted \
	-loj | \
	annotate_regions.py \
		$percent_camo $biotype_camo $coding_camo "camo" \
		> $camo_annotations

## Create GATK bed: bed that will give regions were camo variants will be called
## The GATK bed is the CDS align to regions that are exclusively camo,
## The normal align_to lists the whole genebody element, GATK restricts to just those camo regions
grep -vE "^#" $camo_annotations | \
	awk '$5 == "CDS"' | \
	bedtools intersect \
		-a $alignto_sorted \
		-b -\
	> $gatk_bed
		
# Create bed for regions that are low_mapq but NOT camouflaged
echo "`date` bedtools subtract -a $mapq_merged -b $camo_sorted > $mapq_not_camo"
bedtools subtract \
	-a $mapq_merged \
	-b $camo_sorted \
	> $mapq_not_camo

mkdir -p $RESULT_DIR
cp $mapq_merged $depth_merged $dark_merged \
	$percent_dark $percent_depth $percent_mapq \
	$biotype_dark $biotype_mapq $biotype_depth \
	$dark_annotations $mapq_annotations $depth_annotations \
	$alignto_sorted $realign_sorted $mapq_not_camo \
	$camo_annotations $percent_camo $biotype_camo \
	$coding_dark $coding_mapq $coding_depth $coding_camo \
	$gatk_bed \
	$RESULT_DIR

###################
#### CLEAN UP #####
###################

cd $HOME #move to safe location
rm -rfv $TMP_DIR
echo "---"

echo "`date` DONE"
