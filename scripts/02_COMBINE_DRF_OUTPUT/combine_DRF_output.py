#!/usr/bin/python
# USAGE: python combine_spike_output.py DRF_output_list combined_low_depth combined_low_mapq
#	DRF_output_list: list of DRF outputs to combine
#	combined_low_depth: output file to write the combined low_depth.dark.bed file
#	combined_low_mapq: output file ot write the combined low_mapq.dark.bed file
#
#Must run DRF with --min-depth == -1 and --min-mapq-mass == -1
# This will print out every position into low_mapq.dark.bed to accurately combine output
#
# combined DRF beds will use defualt cutoffs --min-depth = 5 and --min-mapq-mass = 90
import sys
from collections import defaultdict
MIN_DEPTH = 5
MIN_MAPQ_MASS = 90

def main(DRF_bed_list_file, combined_low_depth, combined_low_mapq):
	try:
		beds_to_combine = []
		DRF_bed_list = open(DRF_bed_list_file, 'r')
		for line in DRF_bed_list:
			bed_file_name = line.strip()
			beds_to_combine.append(open(bed_file_name, 'r'))

		combined_low_depth = open(combined_low_depth, 'w')
		combined_low_mapq = open(combined_low_mapq, 'w')
		for line in beds_to_combine[0]:
			if line.startswith("chrom\tstart\tend"):
				for bed in beds_to_combine[1:]:
					line = bed.readline()
					if not line.startswith("chrom\tstart\tend"):
						raise ValueError("Input bed files do noy have headers in same order")
				continue

			toks = line.strip().split('\t')
			chrom = toks[0]
			start = int(toks[1])
			if start % 1000000 == 0:
				sys.stderr.write("Analyzed %d positions on contig %s; continuing\n" % (start, chrom))
			nMapQBelowThresh_total = int(toks[3])
			depth_total = int(toks[4])
			for bed in beds_to_combine[1:]:
				line = bed.readline()
				toks = line.strip().split('\t')
				bed_chrom = toks[0]
				bed_start = int(toks[1])
				if bed_chrom != chrom or bed_start != start:
					raise ValueError("Input bed files are not in the same order")

				nMapQBelowThresh_total += int(toks[3])
				depth_total += int(toks[4])

			avgMapQBelowThresh = float(nMapQBelowThresh_total) / len(beds_to_combine)
			avgDepth = float(depth_total) / len(beds_to_combine)
			avgPercMapQBelowThresh = avgMapQBelowThresh / avgDepth * 100 if avgDepth > 0 else -1.0

			out_line = "%s\t%d\t%d\t%f\t%f\t%f\n" % (chrom, start, start + 1, avgMapQBelowThresh, avgDepth,
					avgPercMapQBelowThresh)
			
			if avgDepth <= MIN_DEPTH:
				combined_low_depth.write(out_line)
			elif avgPercMapQBelowThresh >= MIN_MAPQ_MASS:
				combined_low_mapq.write(out_line)

	except ValueError as e:
		print(e)
			
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])
