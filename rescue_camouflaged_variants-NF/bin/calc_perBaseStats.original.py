#!/usr/bin/env python3

import sys
import mgzip




chrlist = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
chrlist_noChr = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",  "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"]

def main(inputGenomeBed, outputStatFile):

    try:
        cumulativeStats = {} # this is formatted to be {chr: [Cumulative_coverage_sum, Number_of_bases]}
        with mgzip.open(inputGenomeBed, "rt", thread=4) as genomeBed:
            for line in genomeBed:
                line = line.strip("\n").split("\t")
                if line[0] in chrlist or line[0] in chrlist_noChr:
                    if line[0] in cumulativeStats:
                        
                        cumulativeStats[line[0]][0] += float(line[4])
                        cumulativeStats[line[0]][1] += float(line[2])-float(line[1])
                            
                    else:
                        cumulativeStats[line[0]] = [float(line[4]), float(line[2])-float(line[1])]

                    
                    if "AllChrs" in cumulativeStats:
                        cumulativeStats["AllChrs"][0] += float(line[4])
                        cumulativeStats["AllChrs"][1] += float(line[2])-float(line[1])
                         
                    else:
                        cumulativeStats["AllChrs"] = [float(line[4]), float(line[2])-float(line[1])]

        ## Write dict to output
        header="chr\tcumSumDepth\tNumBasePairs\tAvgDepth\n"
        with open(outputStatFile, 'wt') as output:
            output.write(header)
            for chrom in cumulativeStats:
                avgDepth = cumulativeStats[chrom][0]/float(cumulativeStats[chrom][1])
                outline = "%s\t%d\t%d\t%f\n" % (chrom, cumulativeStats[chrom][0], cumulativeStats[chrom][1], avgDepth)
                output.write(outline)

    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
