
# Enable bash debugging to log all commands
set -x

echo "Job running on: `hostname`"

CAMO_ANNOTATION=$1
INCLUDE_REGIONS=$2
MASK_BED=$3
REF=$4
OUT_FILE=$5
NUM_THREADS=$6

query="./query/mapq_genes.query.fa"
blat_result="./blat_result/blat.results.psl"
blat_log="./blat_log/tmp.blat.log"
blat_bed=${blat_result//psl/bed}
false_positives="false_positives.txt"

################################################
# Create BLAT queries from camouflaged regions #
################################################

mkdir -p "query"
if [[ "${INCLUDE_REGIONS,,}" == "cds" ]]; then

	# Only awk out CDS regions (i.e., only include CDS)
	if ! grep -vE "^#" $CAMO_ANNOTATION | \
		awk '$5 == "CDS"' | \
		bedtools intersect \
			-a - \
			-b $MASK_BED \
			-wb | \
			awk '$NF <= 5' | \
			bedtools getfasta \
				-fi $REF \
				-bed - \
				-name+ \
				-fo $query; then
		echo "`date` bedtools intersect failed for $CAMO_ANNOTATION and $MASK_BED"
		exit 1
	fi
else

	# Do NOT awk out CDS regions (i.e., use all)
	if ! grep -vE "^#" $CAMO_ANNOTATION | \
		bedtools intersect \
			-a - \
			-b $MASK_BED \
			-wb | \
			awk '$NF <= 5' | \
			bedtools getfasta \
				-fi $REF \
				-bed - \
				-name+ \
				-fo $query; then
		echo "`date` bedtools intersect failed for $CAMO_ANNOTATION and $MASK_BED"
		exit 1
	fi
fi

nLines=$(wc -l $query | awk '{print $1}')
stepSize=$(($nLines / $NUM_THREADS))
if [[ $(($stepSize % 2 )) == 1 ]]
then
	    stepSize=$(($stepSize + 1))
fi

rm -rf "./blat_result"
mkdir -p "./blat_result"
mkdir -p "./blat_log"
for i in $(seq 1 $stepSize $nLines)
do
	sed "$((i)),$((i + $stepSize - 1))!d" $query > ${query}.${i}
	if ! blat $REF ${query}.${i} \
			-t=dna \
			-q=dna \
			-minIdentity=95 \
			-noHead \
			${blat_result}.${i} \
			&> ${blat_log}.${i} &; then
	echo "`date` blat failed for $REF and ${query}.${i}"
	exit 1
fi
done
wait

echo "`date` blatting complete: combining blat output"
cat ${blat_result}.* > $blat_result

echo "`date` scoring blat output"
if ! awk -f score_blat_output.awk \
	$blat_result \
	> $blat_bed; then
	echo "`date` score_blat_output.awk failed for $blat_result."
	exit 1
fi

echo "`date` extracting false positives"
if ! python extract_false_positives.py \
	$blat_bed \
	$REF | \
	bedtools sort -i - -g ${REF}.fai \
	> $false_positives; then
	echo "`date` extract_false_positives.py failed for $blat_bed"
	exit 1
fi

cp $false_positives $OUT_FILE
