#!/usr/bin/env python3
import sys


valid_chroms = set(
                    list(map(str, range(1,23))) + \
                        ['X', 'Y', 'M', 'MT'] + \
                        ['chr' + chrom for chrom in map(str, range(1,23))] + \
                        ['chrX', 'chrY', 'chrM', 'chrMT']
                   )



for line in sys.stdin:
	toks = line.strip().split('\t')
	chrom = toks[0]
	if chrom == "chrom": continue
	elif chrom.startswith('chr'):
		chrom = chrom[3:]

	if chrom not in valid_chroms:
		#as soon as we hit an unassembled contig, quit
		break 
	sys.stdout.write(line)
