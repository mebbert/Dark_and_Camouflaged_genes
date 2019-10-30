#!/usr/bin/python
import sys
from collections import defaultdict

def extractAnnos(annos_line):
	annos_toks = annos_line.split(';')
	annos = {}
	for annos_tok in annos_toks:
		toks = annos_tok.split('=')
		annos[toks[0]] = toks[1]
	return annos

def getClusterRepresentative(cluster_list):
	max_gene = None
	max_length = 0
	for region in cluster_list:
		length = region[1] - region[0]
		if region[2] == "gene" and length > max_length:
			max_gene = region
			max_length = length
	
	if max_gene: return max_gene
	max_region = None
	for region in cluster_list:
		length = region[1] - region[0]
		if length > max_length:
			max_region = region
			max_length = length
	return max_region
		
def selectGenesToKeep():
	clusters = defaultdict(list)
	for line in sys.stdin:
		toks = line.strip().split("\t")
		region_type = toks[2]
		start = int(toks[3])
		end = int(toks[4])
		annos = extractAnnos(toks[8])
		ID = annos['ID']
		value = (start, end, region_type, ID)
		clusters[toks[-1]].append(value)

	genesToKeep = set()
	for key in clusters:
		clusterRep = getClusterRepresentative(clusters[key])
		genesToKeep.add(clusterRep[-1])
		
	return genesToKeep

def main(gff3_file):
	genesToKeep = selectGenesToKeep()
	gff3 = open(gff3_file, 'r')
	keep_gene = False
	keep_transcript = False
	gene_biotype = ''
	for line in gff3:
		if line.startswith('#'): continue
		toks = line.strip().split('\t')
		chrom = toks[0]
		region_type = toks[2]
		start = int(toks[3])
		end = int(toks[4])
		score = toks[5]
		strand = toks[6]

		annos = toks[8]
		annos_dict = extractAnnos(annos)
		ID = ''
		biotype = ''
		if "ID" in annos_dict: 
			ID = annos_dict["ID"]
		if "biotype" in annos_dict: 
			biotype = annos_dict["biotype"]

		if ID.startswith("gene:"):
			keep_gene = False
			if biotype == "antisense":
				keep_gene = False
			elif ID in genesToKeep:
				keep_gene = True
				keep_transcript = True
				gene_biotype = biotype

		elif keep_gene and ID.startswith("transcript:"):
			keep_transcript = True
			if biotype == "retained_intron" or (gene_biotype == "protein_coding" and biotype != "protein_coding"):
				keep_transcript = False

		if keep_gene and keep_transcript:
			sys.stdout.write(line)

	gff3.close()

if __name__ == "__main__":
	main(sys.argv[1])
