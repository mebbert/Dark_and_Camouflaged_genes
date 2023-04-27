#!/usr/bin/env python3
import itertools
import copy
import sys
from collections import defaultdict

def extractAnnos(annos):
	anno_dict = {}
	anno_toks = annos.split(';')
	for anno in anno_toks:
		toks = anno.split('=')
		anno_dict[toks[0]] = toks[1]

	return anno_dict

def main(percent_camo_file, biotype_camo_file, coding_camo_file, label):
	lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: (0,0))))
	coding_lengths = defaultdict(lambda:defaultdict(lambda: (0,0)))
	biotype_lengths = defaultdict(lambda:defaultdict(lambda: (0,0)))
	gene_names = {}
	biotypes = {}
	biotype_bins = ["total", "protein coding", "pseudogene", "rRNA", "snRNA",
		"miRNA",  "lincRNA", "other"]
	sys.stdout.write("# This table shows the intersection of %s regions by a gene annotation gff3, showing where %s regions fall within gene body elements\n" % (label, label))
	sys.stdout.write("# chrom : Chromosome of %s region\n" % label)
	sys.stdout.write("# region_start : start position of %s region within gene-body element\n" % label)
	sys.stdout.write("# region_end : end position of %s region within gene-body element\n" % label)
	sys.stdout.write("# gene_body_id : ID of gene body element that contains this %s region\n" % label)
	sys.stdout.write("# region_type : genebody element type (e.g. exon, intron, UTR, etc.)\n")
	sys.stdout.write("# biotype : GENCODE biotype of gene that contains this %s reigon\n" % label)
	sys.stdout.write("# gene_body_chrom : Chromosome of the genebody element\n")
	sys.stdout.write("# gene_body_start : start position of the genebody element\n")
	sys.stdout.write("# gene_body_end : end position of the genebody element\n")
	sys.stdout.write("#chrom\tstart\tend\tgene_body_id\tregion_type\tbiotype\tgene_body_chrom\tgene_body_start\tgene_body_end\n")
	for line in sys.stdin:
		toks = line.strip().split('\t')

		chrom = toks[0]
		region_start = int(toks[1])
		region_end = int(toks[2])
		region_length = region_end - region_start

		region_type = toks[3]
		region_bin = region_type
		if "UTR" in region_type: region_bin = "UTR"
		elif region_type not in ["exon", "intron", "CDS"]: region_bin = "gene"

		annos = extractAnnos(toks[6])
		biotype = annos['biotype']
		region_id = annos['ID']
		region_name = annos['Name']
		
		gene_id = "_".join(region_id.split('_')[0:-2])
		if region_type == "gene":
			gene_id = region_id
		if gene_id.startswith("gene"): gene_id = gene_id[5:]
		gene_name = "_".join(region_name.split('_')[0:-2])
		gene_names[gene_id] = gene_name
		
		biotype_bin = biotype
		if "pseudo" in biotype: biotype_bin = "pseudogene"
		elif biotype in ["snRNA", "snoRNA"]: biotype_bin = "snRNA"
		elif biotype in ["lncRNA", "lincRNA"]: biotype_bin = "lincRNA"
		elif biotype not in biotype_bins: biotype_bin = "other"
		if biotype == "protein_coding": biotype_bin = "protein coding"
		biotypes[gene_id] = biotype_bin

		overlap = 0
		if toks[7] != ".":
			camo_start = int(toks[8])
			camo_end = int(toks[9])

			intersect_start = max(camo_start, region_start)
			intersect_end = min(camo_end, region_end)
			overlap = intersect_end - intersect_start
			if region_bin != "gene" and overlap > 22:
				value = (chrom, intersect_start, intersect_end, region_name, region_type, biotype, chrom, region_start, region_end)
				sys.stdout.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n" % value)

		lengths[region_bin][gene_id][region_id] = (lengths[region_bin][gene_id][region_id][0] + overlap, region_length)
		if region_bin in ["CDS", "UTR"]:
			lengths["exon"][gene_id][region_id] = (lengths["exon"][gene_id][region_id][0] + overlap, region_length)
		elif region_bin == "gene":
			biotype_lengths["total"][region_id] = (biotype_lengths["total"][region_id][0] + overlap, region_length)
			biotype_lengths[biotype_bin][region_id] = (biotype_lengths[biotype_bin][region_id][0] + overlap, region_length)

		if biotype_bin == "protein coding":
			coding_lengths[region_bin][region_id] = (coding_lengths[region_bin][region_id][0] + overlap, region_length)
			if region_bin in ["CDS", "UTR"]:
				coding_lengths["exon"][region_id] = (coding_lengths["exon"][region_id][0] + overlap, region_length)

	biotype_camo = open(biotype_camo_file, 'w')
	biotype_camo.write("# This table quanitifies the number of %s dark bases stratified by GENCODE biotype\n" % label)
	biotype_camo.write("# biotype : GENCODE biotype\n")
	biotype_camo.write("# num_dark_bases : number of %s bases in biotype\n" % label)
	biotype_camo.write("# total_biotype_bases : total count of all bases in biotype\n")
	biotype_camo.write("# perc_region : percentage of bases that are %s in this biotype (num_dark_bases / total_biotype_bases) \n" % label)
	biotype_camo.write("biotype\tnum_dark_bases\ttotal_biotype_bases\tperc_dark\n")
	for biotype in biotype_bins:
		total_camo = 0
		total_length = 0
		genes = biotype_lengths[biotype]
		for gene in genes:
			total_camo += genes[gene][0]
			total_length += genes[gene][1]
		if total_length == 0: continue
		biotype_camo.write("%s\t%d\t%d\t%f\n" % (biotype, total_camo, total_length, float(total_camo)/total_length*100))

	percent_camo = open(percent_camo_file, 'w')
	percent_camo.write("# This table shows what percentage of genes are %s, stratefied by genebody element\n" % label)
	percent_camo.write("# gene_name : Ensembl gene name\n")
	percent_camo.write("# biotype : GENCODE biotype of gene\n")
	percent_camo.write("# biotype : percent of CDS bases that are %s within gene (non protein-coding genes will have value of -1)\n" % label)
	percent_camo.write("# perc_UTR : percent of UTR bases that are %s within gene (non protein-coding genes will have value of -1(\n" % label)
	percent_camo.write("# perc_exon : percent of exonic bases that are %s\n" % label)
	percent_camo.write("# perc_intron : percent of intronic bases that are %s\n" % label)
	percent_camo.write("# perc_total : percent of total gene that is %s\n" % label)
	percent_camo.write("gene_name\tbiotype\tperc_CDS\tperc_UTR\tperc_exon\tperc_intron\tperc_total\n")
	for gene_id in gene_names:
		results = defaultdict(lambda: -1)
		for region_bin in lengths:
			region_lengths = lengths[region_bin][gene_id]
			total_camo = 0
			total_length = 0
			for region in region_lengths:
				total_camo += region_lengths[region][0]
				total_length += region_lengths[region][1]
			if total_length == 0: continue
			results[region_bin] = float(total_camo/total_length)*100
		name = gene_names[gene_id]
		biotype = biotypes[gene_id]
		if results["gene"] == 0.0: continue
		percent_camo.write("%s\t%s\t%f\t%f\t%f\t%f\t%f\n" % (name, biotype, results["CDS"], results["UTR"], results["exon"], results["intron"], results["gene"]))

	coding_camo = open(coding_camo_file, 'w')
	coding_camo.write("# This table quantifies where dark regions fall within protein-coding genes, stratifying the percentage of %s bases by gene-body element\n" % label)
	coding_camo.write("# region_type : gene-body element ( e.g. Exon, Intron, CDS, UTR )\n")
	coding_camo.write("# count : total number of %s bases within this gene-body element\n" % label)
	coding_camo.write("# total : total count of bases in across all gene-body elements of this type\n")
	coding_camo.write("region_type\tcount\ttotal\tperc_region\n")
	for region_type in coding_lengths:
		total_camo = 0
		total_length = 0
		regions = coding_lengths[region_type]
		for region_id in regions:
			total_camo += regions[region_id][0]
			total_length += regions[region_id][1]
		if total_length == 0: continue
		coding_camo.write("%s\t%d\t%d\t%f\n" % (region_type, total_camo, total_length, float(total_camo)/total_length*100))
	coding_camo.close()

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
