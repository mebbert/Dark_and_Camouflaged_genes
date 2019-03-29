#!/usr/bin/python
import sys
from collections import defaultdict

genes = defaultdict(int)
for line in sys.stdin:
	toks = line.strip().split(';')
	name_tok = toks[1]
	gene_region = name_tok.split('=')[1]
	gene = gene_region.split('_')[0]
	genes[gene] += 1

print("gene\tvariant_count")
print("total_genes\t%d" % len(genes))
total = 0
for gene,count in genes.items():
	print("%s\t%d" % (gene, count))
	total += count
print("total\t%d" % total)

	
