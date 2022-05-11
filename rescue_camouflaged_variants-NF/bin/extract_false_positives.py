#!/usr/bin/python
from Bio import SeqIO, pairwise2
import sys

def getSeqFromRef(RefDict, pos):
    #Blat results are 1-based, pyhton Seq object is 0-based
    if pos[-1] == "-": 
            #if coordinates come from reverse strand, return reverse complement
            return RefDict[pos[0]][(pos[1] - 1):pos[2]].reverse_complement().seq
    else:
            return RefDict[pos[0]][(pos[1] - 1):pos[2]].seq


def main(blat_bed, ref):
    RefDict = SeqIO.index(ref, "fasta")
    blat_bed = open(blat_bed, 'r')
    seen_variants = set()
    for line in blat_bed:
        toks = line.strip().split()
        target_pos = (toks[0], int(toks[1]), int(toks[2]), toks[-1])
        query_tok = toks[3].split(":")
        start, end = query_tok[-1].split("-")
        query_pos = (query_tok[-2], int(start), int(end), "+") #Query always on forward strand
        target = getSeqFromRef(RefDict, target_pos)
        query = getSeqFromRef(RefDict, query_pos)
        for a in pairwise2.align.localms(target, query, 1, -3, -5, -2):
            #sys.stderr.write(pairwise2.format_alignment(*a))
            aligned_target = a[0]
            aligned_query = a[1]
            i = 0
            end = len(aligned_target) - 1

            while aligned_target[i]  == "-" or aligned_query[i] == "-": i += 1
            while aligned_target[end] == "-" or aligned_query[end] == "-": end -= 1

            while i <= end:
                if aligned_target[i] == "-":
                    start_del = i - 1 # in VCFs indel positions 1 base before actual indel
                    while aligned_target[i] == "-":
                        i += 1
                    var_start = query_pos[1] + start_del
                    var = "%s\t%d\t%d\t%s\t%s" % (query_pos[0], var_start, var_start + 1,
                            aligned_query[start_del:i], aligned_target[start_del])
                    if var not in seen_variants:
                        seen_variants.add(var)
                        print(var)
                    break

                if aligned_query[i] == "-":
                    start_ins = i - 1 # in VCF indel positions are 1 base before actual indel
                    while aligned_query[i] == "-":
                        i += 1
                    var_start = query_pos[1] + start_ins
                    var = "%s\t%d\t%d\t%s\t%s" % (query_pos[0], var_start, var_start + 1,
                            aligned_query[start_ins], aligned_target[start_ins:i])
                    if var not in seen_variants:
                        seen_variants.add(var)
                        print(var)
                        sys.stderr.write(var + "\n")
                    break

                if aligned_query[i] != aligned_target[i]:
                    var_start = query_pos[1] + i
                    var = "%s\t%d\t%d\t%s\t%s" % (query_pos[0], var_start, var_start + 1, 
                                    aligned_query[i], aligned_target[i])
                    if var not in seen_variants:
                        seen_variants.add(var)
                        print(var)
                        sys.stderr.write(var + "\n")

                i += 1

            break
                                    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])


