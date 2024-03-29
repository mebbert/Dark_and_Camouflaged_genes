#!/usr/bin/python
import sys
from collections import defaultdict

def extractAnnos(annos):
	annos_dict = {}
	annos_toks = annos.split(';')
	for annos_tok in annos_toks:
		toks = annos_tok.split("=")
		annos_dict[toks[0]] = toks[1]

	return annos_dict

valid_regions = set(list(map(str, range(1,23))) + ['X', 'Y', 'MT'])
def printLine(region):
	if region == None:
		return
	elif region[0] in valid_regions:
		print("%s\t%d\t%d\t%s\t%s\t%s\tID=%s;Name=%s;biotype=%s" % region)
		
def main():
	current_gene = ""
	previous_region = None
	exon_count = defaultdict(int)
	intron_count = defaultdict(int)
	UTR_count = defaultdict(int)
	gene_beginning = False
	for line in sys.stdin:
		toks = line.strip().split('\t')
		chrom = toks[0]
		gene_start = int(toks[1])
		gene_end = int(toks[2])
		gene_type = toks[3]
		gene_score = toks[4]
		gene_strand = toks[5]
		annos = extractAnnos(toks[6])
		gene_ID = annos["ID"]
		gene_name = annos["Name"]
		biotype = annos["biotype"]

		if toks[7] == '.': 
			continue

		if current_gene != gene_ID:
			printLine(previous_region)
			current_gene = gene_ID
			gene = (chrom, gene_start, gene_end, gene_type, gene_score, gene_strand, gene_ID, gene_name, biotype)
			printLine(gene)
			previous_region = None
			gene_beginning = True

		region_start = int(toks[8])
		region_start = max(gene_start, region_start)
		region_end = int(toks[9])
		region_end = min(gene_end, region_end)
		region_type = toks[10]
		if ',' in region_type: region_type = "UTR"
		if biotype == "protein_coding" and region_type == "exon": region_type = "CDS"
		region_score = toks[11]
		region_strand = toks[12]
		region_ID = ""
		region_Name = ""
		if "UTR" in region_type:
			UTR_count[gene_ID] += 1
			region_ID = gene_ID[5:] + "_UTR_%d" % UTR_count[gene_ID]
			region_Name = gene_name + "_UTR_%d" % UTR_count[gene_ID]
		else:
			exon_count[gene_ID] += 1
			region_ID = gene_ID[5:] + "_%d" % exon_count[gene_ID]
			region_Name = gene_name + "_%d" % exon_count[gene_ID]

		region = (chrom, region_start, region_end, region_type, region_score, region_strand, region_ID, region_Name, biotype)

		if gene_beginning and not previous_region:
			previous_region = region
		elif gene_beginning and "UTR" in region_type and previous_region[3] == "CDS":
			gene_beginning = False
			previous_region = (chrom, previous_region[1], region[2], region_type, region_score, region_strand, region_ID, region_Name, biotype) 
			exon_count[gene_ID] -= 1
		else:
			gene_beginning = False
			printLine(previous_region)

			#check to see if there is gaps between exons (if there is fill with intron)
			if region[1] - previous_region[2] - 1 > 0:
				intron_start = previous_region[2] + 1
				intron_end = region[1] - 1
				intron_count[gene_ID] += 1
				intron_ID = gene_ID[5:] + "_intron_%d" % intron_count[gene_ID]
				intron_Name = gene_name + "_intron_%d" % intron_count[gene_ID]
				intron = (chrom, intron_start, intron_end, "intron", region_score, region_strand, intron_ID, intron_Name, biotype)
				printLine(intron)

			previous_region = region
	printLine(previous_region)

if __name__ == "__main__":
	main()
