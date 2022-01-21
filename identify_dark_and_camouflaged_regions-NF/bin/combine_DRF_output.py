#!/usr/bin/env python
# USAGE: python combine_DRF_output.py low_mapq_file_list combined_low_depth_out_file
# combined_low_mapq_out_file
#	low_mapq_file_list: list of DRF (low_mapQ) outputs to combine
#	combined_low_depth_out_file: output file to write the combined low_depth.dark.bed file
#	combined_low_mapq_out_file: output file ot write the combined low_mapq.dark.bed file
#
# For individual samples, DRF must be run with --min-depth == -1 and --min-mapq-mass == -1.
# This will print out every position into low_mapq.dark.bed to accurately combine output
#
# For this script, combined low_mapq .beds will use default cutoffs --min-depth = 5 and --min-mapq-mass = 90

import sys
from collections import defaultdict

# The minimum depth for a given nucleotide (and run of nucleotides)
# to NOT be considered 'dark-by-depth'
MIN_DEPTH = 5

# The minimum percentage (i.e., mass) of reads aligned to a given nucleotide
# to be considered 'dark-by-mapping quality
MIN_MAPQ_MASS = 90

def main(low_mapq_file_list, combined_low_depth_out_file, combined_low_mapq_out_file):
    try:
        beds_to_combine = []
        low_mapq_bed_list = open(low_mapq_file_list, 'r')
        print("opened the low_mapq_bed file list")
        for line in low_mapq_bed_list:
            bed_file_name = line.strip()
            beds_to_combine.append(open(bed_file_name, 'r'))

        combined_low_depth = open(combined_low_depth_out_file, 'w')
        combined_low_mapq = open(combined_low_mapq_out_file, 'w')
        print(beds_to_combine)

        for line in beds_to_combine[0]:

            # Verify all input bed files have the same header line
            if line.startswith("chrom\tstart\tend"):
                for bed in beds_to_combine[1:]:
                    line = bed.readline()
                    if not line.startswith("chrom\tstart\tend"):
                        raise ValueError("Input bed files do not have headers in same order")
                continue

            toks = line.strip().split('\t')
            chrom = toks[0]
            start = int(toks[1])
            if start % 1000000 == 0:
                sys.stderr.write("Analyzed %d positions on contig %s; continuing\n" % (start, chrom))

            # The number of reads with a MAPQ less than or equal to the defined '--mapq_threshold'
            # when DRF was run. The recommended GATK threshold is < 10, and is what we use as the 
            # default for our pipeline (except DRF uses '<=', so we use 9 as the input).
            nMapQBelowThresh_total = int(toks[3])
            depth_total = int(toks[4])
            for bed in beds_to_combine[1:]:
                line = bed.readline()
                toks = line.strip().split('\t')
                bed_chrom = toks[0]
                bed_start = int(toks[1])

                # All input .bed files should have all positions for the reference
                # genome, so we should never be out of sync. Fail if we are.
                if bed_chrom != chrom or bed_start != start:
                    raise ValueError("Input bed files are not in the same order")

                nMapQBelowThresh_total += int(toks[3])
                depth_total += int(toks[4])

            avgMapQBelowThresh = float(nMapQBelowThresh_total) / len(beds_to_combine)
            avgDepth = float(depth_total) / len(beds_to_combine)
            avgPercMapQBelowThresh = avgMapQBelowThresh / avgDepth * 100 if avgDepth > 0 else -1.0

            out_line = "%s\t%d\t%d\t%f\t%f\t%f\n" % (chrom, start, start + 1, avgMapQBelowThresh, avgDepth, avgPercMapQBelowThresh)
                    
            if avgDepth <= MIN_DEPTH:
                combined_low_depth.write(out_line)
            elif avgPercMapQBelowThresh >= MIN_MAPQ_MASS:
                combined_low_mapq.write(out_line)

    except ValueError as e:
        print(e)
			
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
