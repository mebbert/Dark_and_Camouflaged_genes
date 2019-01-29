#!/usr/bin/python
import sys

def extractAnnos(annos):
	anno_dict = {}
	anno_toks = annos.split(';')
	for anno in anno_toks:
		toks = anno.split('=')
		anno_dict[toks[0]] = toks[1]

	return anno_dict

def main(camo_bed, hg38_annotation):
	region_dict = {}
	hg38_annotation = open(hg38_annotation, 'r')
	for line in hg38_annotation:
		toks = line.strip().split("\t")
		chrom = toks[0]
		region_start = toks[1]
		region_end = toks[2]
		
		annos = extractAnnos(toks[6])
		region_id = annos['Name']
		region_dict[region_id] = (chrom, region_start, region_end)
	hg38_annotation.close()

	camo_bed = open(camo_bed, 'r')
	for line in camo_bed:
		toks = line.strip().split('\t')
		region_id = toks[3].split(';')[0]
		if region_id not in region_dict:
			continue
		new_pos = region_dict[region_id]
		toks[0] = new_pos[0]
		toks[1] = new_pos[1]
		toks[2] = new_pos[2]
		print('\t'.join(toks))
	camo_bed.close()

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
