#!/usr/bin/env python3
import sys
from collections import defaultdict

def extractAnnos(annos):
    annos_dict = {}
    annos_toks = annos.split(';')
    for annos_tok in annos_toks:
        toks = annos_tok.split("=")
        annos_dict[toks[0]] = toks[1]

    return annos_dict

valid_regions = set(
                    list(map(str, range(1,23))) + \
                        ['X', 'Y', 'M', 'MT'] + \
                        ['chr' + chrom for chrom in map(str, range(1,23))] + \
                        ['chrX', 'chrY', 'chrM', 'chrMT']
                   )

def printLine(region):
    if region == None:
        return
    elif region[0] in valid_regions:
        print("%s\t%d\t%d\t%s\t%s\t%s\tID=%s;Name=%s;biotype=%s" % region)
    #print(region[0])
        
def main():
    current_gene = ""
    previous_region = None
    exon_count = defaultdict(int)
    intron_count = defaultdict(int)
    UTR_count = defaultdict(int)
    gene_beginning = False
    for line in sys.stdin:

        # Collect meta data
        toks = line.strip().split('\t')
        chrom = toks[0]
        gene_start = int(toks[1])
        gene_end = int(toks[2])
        gene_type = toks[3]
        gene_score = toks[4]
        gene_strand = toks[5]
        annos = extractAnnos(toks[6])

        # We expect annos["ID"] to return either "ID=<ENSG_ID>" or "ID=gene:<ENSG_ID>". Strip
        # "gene:" if it exists.
        gene_ID = annos["ID"]
        if gene_ID.startswith('gene:'):
            gene_ID = gene_ID.split(':')[1]

        biotype = annos["biotype"]

        # Somewhere between GRCh38 release 93 and 105, Ensembl dropped 'Name' annotations for some
        # genes. Not sure why. If the 'Name' annotation is missing, just set it to the gene ID.
        try:
            gene_name = annos["Name"]
        except KeyError as e:
            gene_name = gene_ID

        if toks[7] == '.': 
           # print('toks 7 is a dot')
            continue

        # Entering a new gene definition block in the gff file
        if current_gene != gene_ID:
            #print('current_gene != gene_ID')
            printLine(previous_region)
            current_gene = gene_ID
            gene = (chrom, gene_start, gene_end, gene_type, gene_score, gene_strand, gene_ID, gene_name, biotype)
            #print("printing gene")
            #print(gene)
            printLine(gene)
            previous_region = None
            gene_beginning = True

        # Collect region information
        region_start = int(toks[8])
        region_start = max(gene_start, region_start)
        region_end = int(toks[9])
        region_end = min(gene_end, region_end)
        #if region_end == -1: print("-1 region end") 
        region_type = toks[10]
        if ',' in region_type: 
            region_type = "UTR"
        if biotype == "protein_coding" and region_type == "exon": 
            region_type = "CDS"
        region_score = toks[11]
        region_strand = toks[12]
        region_ID = ""
        region_Name = ""

        # Handle UTR regions
        if "UTR" in region_type:
            #print('utr in region_type')
            UTR_count[gene_ID] += 1
            region_ID = gene_ID + "_UTR_%d" % UTR_count[gene_ID]
            region_Name = gene_name + "_UTR_%d" % UTR_count[gene_ID]
        else:
            #print('utr NOT in region_type')
            exon_count[gene_ID] += 1
            region_ID = gene_ID + "_exon_%d" % exon_count[gene_ID]
            region_Name = gene_name + "_exon_%d" % exon_count[gene_ID]

        region = (chrom, region_start, region_end, region_type, region_score, region_strand, region_ID, region_Name, biotype)

        if gene_beginning and not previous_region:
            previous_region = region
        elif gene_beginning and "UTR" in region_type and previous_region[3] == "CDS":
            gene_beginning = False
            previous_region = (chrom, previous_region[1], region[2], region_type, region_score, region_strand, region_ID, region_Name, biotype) 
            exon_count[gene_ID] -= 1
        else:
            gene_beginning = False
            #print("printing previous region")
            printLine(previous_region)

            #check to see if there are gaps between exons (if there are fill with intron)
            if region[1] - previous_region[2] - 1 > 0:
                intron_start = previous_region[2] + 1
                intron_end = region[1] - 1
                intron_count[gene_ID] += 1
                intron_ID = gene_ID + "_intron_%d" % intron_count[gene_ID]
                intron_Name = gene_name + "_intron_%d" % intron_count[gene_ID]
                intron = (chrom, intron_start, intron_end, "intron", region_score, region_strand, intron_ID, intron_Name, biotype)
                #print("printing intron")
                printLine(intron)

            previous_region = region
    #print("previous region")
    printLine(previous_region)

if __name__ == "__main__":
    main()
