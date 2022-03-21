#!/bin/bash

# Enable bash debugging to log all commands
set -x

echo "`date` Running camo_gene_pipeline.sh"

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
			echo "ERROR: unrecognized command line argument: $key; exiting"
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

# if [[ -z $RESULT_DIR ]]
# then
# 	RESULT_DIR=${GVERS}/${SEQUENCER}
# fi

# TMP_DIR="tmp/${SEQUENCER}.${GVERS}/"
# mkdir -p $TMP_DIR

##----Define Variables----------
faidx="${GENOME}.fai"
mapq_mass90="low_mapq-mass90.bed"
mapq_merged="${SEQUENCER}.${GVERS}.low_mapq-merged.bed"
depth_merged="${SEQUENCER}.${GVERS}.low_depth-merged.bed"
dark_merged="${SEQUENCER}.${GVERS}.dark-merged.bed"
mapq_lengths="${SEQUENCER}.${GVERS}.low_mapq.lengths.txt"
dark_loj="dark_loj.txt"
mapq_loj="mapq_loj.txt"
depth_loj="depth_loj.txt"
dark_annotations="${SEQUENCER}.${GVERS}.dark_annotations.txt"
mapq_annotations="${SEQUENCER}.${GVERS}.mapq_annotations.txt"
depth_annotations="${SEQUENCER}.${GVERS}.depth_annotations.txt"
camo_annotations="${SEQUENCER}.${GVERS}.camo_annotations.txt"
percent_dark="${SEQUENCER}.${GVERS}.percent_dark_genes.txt"
percent_mapq="${SEQUENCER}.${GVERS}.percent_low_mapq_genes.txt"
percent_depth="${SEQUENCER}.${GVERS}.percent_low_depth_genes.txt"
percent_camo="${SEQUENCER}.${GVERS}.percent_camo_genes.txt"
biotype_dark="${SEQUENCER}.${GVERS}.percent_dark_biotypes.txt"
biotype_mapq="${SEQUENCER}.${GVERS}.percent_mapq_biotypes.txt"
biotype_depth="${SEQUENCER}.${GVERS}.percent_depth_biotypes.txt"
biotype_camo="${SEQUENCER}.${GVERS}.percent_camo_biotypes.txt"
coding_dark="${SEQUENCER}.${GVERS}.percent_dark_coding_regions.txt"
coding_mapq="${SEQUENCER}.${GVERS}.percent_mapq_coding_regions.txt"
coding_depth="${SEQUENCER}.${GVERS}.percent_depth_coding_regions.txt"
coding_camo="${SEQUENCER}.${GVERS}.percent_camo_coding_regions.txt"
query="query/mapq_genes.query.fa"
blat_result="blat_result/blat.results.psl"
blat_log="blat_log/tmp.blat.log"
blat_bed=${blat_result//psl/bed}
mapped_blat_results="query.results.filtered.mapped.txt"
camo_bed="tmp_camo.bed"
mask_bed="tmp_mask_bed.bed"
extraction="tmp_extraction.bed"
camo_sorted="tmp_camo.sorted.bed"
alignto_sorted="${SEQUENCER}.${GVERS}.camo.mask_bed.sorted.bed"
extraction_sorted="${SEQUENCER}.${GVERS}.camo.extraction.sorted.bed"
gatk_bed_all="${SEQUENCER}.${GVERS}.camo.GATK.all_camo_regions.bed"
gatk_bed_CDS_only="${SEQUENCER}.${GVERS}.camo.GATK.CDS_regions_only.bed"
mapq_gene_bodies="low_mapq_gene_bodies.bed"
mapq_not_camo="${SEQUENCER}.${GVERS}.low_mapq.NOT_camo.bed"

################################
# Function to test exit status #
################################
# function pipe_failed {
# 
# 	exit_status_list=("$@")
# 
# 	for status in "${exit_status_list[@]}"
# 	do
# 		# If exit status is not 0 (success) or 141 (SIGPIPE), throw error. Igore 141 because some
# 		# commands intentionally close the pipe prematurely (e.g., remove_unassembled_contigs.py).
# 		if [[ $status -ne 0 && $status -ne 141 ]]; then
# 			# Return "true" (i.e., 0)
# 			return 0
# 		fi
# 	done
# 
# 	# Return "false" (i.e., 1)
# 	return 1
# }




####################################################################################################
# These first two processes can run simultaneously, but couldn't find a good way to capture failed #
# exit status if we sent them to background with '&'.                                              #
####################################################################################################

###################################################################################
# Merge coordinates for depth bed, removing regions that are less than 20bp long. #
# Calculate the mean, median depth (column 5).                                    #
###################################################################################
# Swallow 141 (SIGPIPE) because remove_unassembled_contigs.py closes pipe abruptly.
# time pigz -dcp 4 $DEPTH_BED | bedtools merge -d 20 -c 5 -o mean,median -i stdin | \
# time bedtools merge -d 20 -c 5 -o mean,median -i $DEPTH_BED | \

# print header line
echo -e "#chrom\tstart\tend\tmean depth\tmedian depth" > ${depth_merged}
time bedtools merge -d 20 -c 5 -o mean,median -i $DEPTH_BED | \
	remove_unassembled_contigs.py | \
	awk '{ if($3 - $2 > 20) print $0}' \
	>> $depth_merged \
	|| if [[ $? -eq 141 ]]; then  # Swallow 141 (SIGPIPE)
			true
		else
			echo "ERROR (`date`): Failed to merge coordinates for $DEPTH_BED. See log for details."
			exit 1
		fi


# Check if the exit status for every cmd in the pipe were successful.
# Cannot simply use 'if ! <cmd>' if piping multiple cmds together. The environment
# variable $PIPESTATUS stores the status for all cmds in the pipe.
TMPSTATUS=("${PIPESTATUS[@]}")
#if pipe_failed "${TMPSTATUS[@]}"; then
#	echo "ERROR (`date`): Failed to merge coordinates for $DEPTH_BED. See log for details."
#	exit 1
#fi
	


##################################################################################
# Merge coordinates for mapq bed, removing regions that are less than 20bp long. #
# Calculate the mean, median depth (column 5).                                   #
##################################################################################
# Swallow 141 (SIGPIPE) because remove_unassembled_contigs.py closes pipe abruptly.
# time bedtools merge -d 20 -c 5 -o mean,median -i $MAPQ_BED | \
# time pigz -dcp 4 $MAPQ_BED | bedtools merge -d 20 -c 5 -o mean,median -i stdin | \

# print header line
echo -e "#chrom\tstart\tend\tmean depth\tmedian depth" > ${mapq_merged}
time bedtools merge -d 20 -c 5 -o mean,median -i $MAPQ_BED | \
	remove_unassembled_contigs.py | \
	awk '{ if($3 - $2 > 20) print $0}' \
	>> $mapq_merged \
	|| if [[ $? -eq 141 ]]; then  # Swallow 141 (SIGPIPE)
			true
		else
			echo "ERROR (`date`): Failed to merge coordinates for $MAPQ_BED. See log for details."
			exit 1
		fi


TMPSTATUS=("${PIPESTATUS[@]}")
#if pipe_failed "${TMPSTATUS[@]}"; then
#	echo "ERROR (`date`): Failed to merge coordinates for $MAPQ_BED. See log for details."
#	exit 1
#fi


# Printing the lengths of each mapq region to a file
echo -e "region_id\tchrom\tlength" >> $mapq_lengths
if ! awk '{n +=1; print n,$1,($3 - $2)}' $mapq_merged >> $mapq_lengths; then
	echo "ERROR: `date` awk failed for $mapq_merged"
	exit 1
fi


########################################
# Merge dark-by-MAPQ and dark-by-depth #
########################################
cat $mapq_merged $depth_merged | \
	bedtools sort -i - -faidx $faidx | \
	bedtools merge -i - \
	> $dark_merged \
	|| {
			echo "ERROR (`date`): Failed to merge dark-by-MAPQ and dark-by-depth. See log for details."
			exit 1
		}

# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to merge dark-by-MAPQ and dark-by-depth. See log for details."
# 	exit 1
# fi


########################################################
## Intersect all dark regions with gene annotation bed #
########################################################
bedtools intersect \
	-a $ANNOTATION \
	-b $dark_merged \
	-loj | \
	annotate_regions.py \
		$percent_dark $biotype_dark $coding_dark "dark" \
		> $dark_annotations \
	|| {
			echo "ERROR (`date`): Failed to intersect all dark regions with gene annotation bed. See log for details."
			exit 1
		}


# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to intersect all dark regions with gene annotation bed. See log for details."
# 	exit 1
# fi




##################################################################
## Intersect the low depth dark regions with gene annotation bed #
##################################################################
bedtools intersect \
	-a $ANNOTATION \
	-b $depth_merged \
	-loj | \
	annotate_regions.py \
		$percent_depth $biotype_depth $coding_depth "dark by depth" \
		> $depth_annotations \
	|| {
			echo "ERROR (`date`): Failed to intersect the low depth dark regions with gene annotation bed. See log for details."
			exit 1
		}

# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to intersect the low depth dark regions with gene annotation bed. See log for details."
# 	exit 1
# fi



############################################################
## Intersect the low mapq regions with gene annotation bed #
############################################################
bedtools intersect \
	-a $ANNOTATION \
	-b $mapq_merged \
	-loj | \
	annotate_regions.py \
		$percent_mapq $biotype_mapq $coding_mapq "dark by MAPQ" \
		> $mapq_annotations \
	|| {
			echo "ERROR (`date`): Failed to intersect the low mapq regions with gene annotation bed. See log for details."
			exit 1
		}


# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to intersect the low mapq regions with gene annotation bed. See log for details."
# 	exit 1
# fi




## Take mapq_file and create a blat query sequence list
#
# NOTE: Somewhere between the bedtools version we originally
# used (v 2.26) and the current version (v 2.30), bedtools changed
# the behavior for `-name` to include the name from the .bed file
# and the coordinates. This feature is useful, but not exactly 
# backwards compatible. We now need to specify `-nameOnly`.
mkdir -p "query"
cat $mapq_annotations | \
	grep -vE "^#" | \
	bedtools getfasta \
		-fi $GENOME \
		-bed - \
		-nameOnly \
		-fo $query \
	|| {
			echo "ERROR (`date`): Failed to create blat query with 'bedtools getfasta'. See log for details."
			exit 1
		}

# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to create blat query with 'bedtools getfasta'. See log for details."
# 	exit 1
# fi



# Run BLAT
nLines=$(wc -l $query | awk '{print $1}')
stepSize=$(($nLines / $NUM_THREADS))
if [[ $(($stepSize % 2 )) == 1 ]]
then
	stepSize=$(($stepSize + 1))
fi

rm -rf "./blat_result"
mkdir -p "./blat_result"
mkdir -p "./blat_log"
pid_array=()
for i in $(seq 1 $stepSize $nLines)
do
	sed "$((i)),$((i + $stepSize - 1))!d" $query > ${query}.${i}
	blat $GENOME ${query}.${i} \
			-t=dna \
			-q=dna \
			-minIdentity=95 \
			-noHead \
			${blat_result}.${i} \
			&> ${blat_log}.${i} &
	
	# Save PID for each submitted background process. "$!" always
	# holds the PID for most recently submitted process.
	pid_array+=($!)
	# echo "\$!: $!"

done
# wait


# Wait until all submitted processes are completed.
keep_waiting=true
while $keep_waiting
do
	keep_waiting=false
	for pid in "${pid_array[@]}"
	do
		if ps -p $pid > /dev/null; then

			# At least this process is still running, so keep waiting
			keep_waiting=true
			echo "Waiting on blat (PID: $pid)..."
		fi
	done

	# sleep for 10 seconds
	sleep 10
done

# Check exit status for each PID and exit if any failed
had_failure=false
for pid in "${pid_array[@]}"
do
	if ! wait $pid; then
		$had_failure=true
		echo "ERROR (`date`): blat failed (PID: $pid)"
	fi
done

if $had_failure; then
	exit 1
fi

## Combine the blat results together and format them
## into a bed file where score shows the sequence identity,
## removing blat hits that are less than 98% sequence identity
# TODO: Should we be including unplaced contigs (e.g., hs38d1)? Probably, 
# but CHM13 incorporates all of them, so will probably just move forward.
echo "`date` blatting complete: combining blat output"
if ! cat ${blat_result}.* > $blat_result; then
	echo "ERROR (`date`): cat failed for ${blat_result}.*"
	exit 1
fi

if ! score_blat_output.awk \
	$blat_result \
	> $blat_bed; then
	echo "ERROR: `date` score_blat_output.awk failed for ${blat_result}"
	exit 1
fi

if ! awk -v "OFS=\t" '$1 !~ "^#" { print $7,$8,$9,$4 }' \
	$mapq_annotations \
	> $mapq_gene_bodies; then
	echo "ERROR: `date` score_blat_output.awk failed for ${blat_result}"
	exit 1
fi


###################################################################
## Intersect the blat output to original list of low_mapq regions #
###################################################################
bedtools intersect \
	-a $blat_bed \
	-b $mapq_gene_bodies \
	-loj | \
	awk '$6 != "."' \
	> $mapped_blat_results \
	|| {
			echo "ERROR (`date`): Failed to Intersect the blat output to original list of low_mapq regions. See log for details."
			exit 1
		}

# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to Intersect the blat output to original list of low_mapq regions. See log for details."
# 	exit 1
# fi



## From the scored intersected blat output, calculates the maps
## of which camo regions align to which other ones, and creates camo sets
# lists all the regions in a camo set in the extraction file and selects the one region form that set that we 
## will use to align to (written in the mask_bed file, all other regions in set will be masked)
if ! extract_camo_regions.py \
	$mapq_annotations \
	$mapped_blat_results \
	$extraction \
	$mask_bed \
	$camo_bed; then
	echo "ERROR: `date` extract_camo_regions.py failed."
	exit 1
fi

echo "`date` sorting camo bed files"

#extract_camo_regions prints out whole camo gene bodies,
#also same gene bodies may be printed twice, so first merge then
#intersect with merged_mapq to get just the Camo Region boundaries
bedtools sort -i $camo_bed -faidx $faidx | \
	bedtools merge -i - | \
	bedtools intersect -a - -b $mapq_merged > $camo_sorted \
	|| {
			echo "ERROR (`date`): bedtools sort/merge/intersect failed for $camo_bed. See log for details."
			exit 1
		}

# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): bedtools sort/merge/intersect failed for $camo_bed. See log for details."
# 	exit 1
# fi



if ! bedtools sort -i $mask_bed -faidx $faidx > $alignto_sorted; then
	echo "ERROR: `date` bedtools sort failed for $mask_bed."
	exit 1
fi
if ! bedtools sort -i $extraction  -faidx $faidx > $extraction_sorted; then
	echo "ERROR: `date` bedtools sort failed for $extraction."
	exit 1
fi

####################################################################################
## Create Camo Annotation table (intersecting camo regions to gene annotation bed) #
####################################################################################
bedtools intersect \
	-a $ANNOTATION \
	-b $camo_sorted \
	-loj | \
	annotate_regions.py \
		$percent_camo $biotype_camo $coding_camo "camo" \
		> $camo_annotations \
	|| {
			echo "ERROR (`date`): Failed to create Camo Annotation table. See log for details."
			exit 1
		}

# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to create Camo Annotation table. See log for details."
# 	exit 1
# fi

## Create GATK bed: bed that will give regions were camo variants will be called
## The GATK bed is the CDS mask_bed regions that are exclusively camo,
## The normal mask_bed lists the whole genebody element, GATK restricts to just those camo regions
grep -vE "^#" $camo_annotations | \
	awk '$5 == "CDS"' | \
	bedtools intersect \
		-a $alignto_sorted \
		-b -\
	> $gatk_bed_CDS_only \
	|| {
			echo "ERROR (`date`): Failed to create GATK bed. See log for details."
			exit 1
		}


## Create GATK bed for all camo regions: bed that will give regions were camo variants will be called
## The GATK bed is for all mask_bed camo regions (not just CDS) that are exclusively camo,
## The normal mask_bed lists the whole genebody element, GATK restricts to just those camo regions
grep -vE "^#" $camo_annotations | \
	bedtools intersect \
		-a $alignto_sorted \
		-b -\
	> $gatk_bed_all \
	|| {
			echo "ERROR (`date`): Failed to create GATK bed. See log for details."
			exit 1
		}



# TMPSTATUS=("${PIPESTATUS[@]}")
# if pipe_failed "${TMPSTATUS[@]}"; then
# 	echo "ERROR (`date`): Failed to create GATK bed. See log for details."
# 	exit 1
# fi



################################################################
# Create bed for regions that are low_mapq but NOT camouflaged #
################################################################
if ! bedtools subtract \
	-a $mapq_merged \
	-b $camo_sorted \
	> $mapq_not_camo; then
	echo "ERROR (`date`): Failed to create bed for regions that are low_mapq but NOT camouflaged. See log for details."
	exit 1
fi

# mkdir -p $RESULT_DIR
# cp $mapq_merged $depth_merged $dark_merged \
# 	$percent_dark $percent_depth $percent_mapq \
# 	$biotype_dark $biotype_mapq $biotype_depth \
# 	$dark_annotations $mapq_annotations $depth_annotations \
# 	$alignto_sorted $extraction_sorted $mapq_not_camo \
# 	$camo_annotations $percent_camo $biotype_camo \
# 	$coding_dark $coding_mapq $coding_depth $coding_camo \
# 	$gatk_bed \
# 	$RESULT_DIR

###################
#### CLEAN UP #####
###################

# cd $HOME #move to safe location
# rm -rfv $TMP_DIR
echo "---"

echo "`date` DONE"

