#!/usr/bin/env python3
import sys
import os
import re
from collections import defaultdict

class Graph:
    def __init__(self):
        self.nodes = defaultdict(set)

    def addEdge(self, a, b):
        self.nodes[a].add(b)

#		if b in self.nodes:
#			self.nodes[b].add(a)
#		else:
#			neighbors = set()
#			neighbors.add(a)
#			self.nodes[b] = neighbors
#
    def DFS(self):
        visited = defaultdict(int)
        DFSmodules = []

        stack = []
        for node in self.nodes:
            if visited[node] == 0:
                DFSmodules.append(set())
                stack.append(node)
                
                while stack:
                    current = stack.pop()
                    DFSmodules[-1].add(current)
                    visited[current] = 1

                    for v in self.nodes[current]:
                        if visited[v] == 0:
                            stack.append(v)

        return DFSmodules

    def printGraph(module):
        print("strict digraph {")
        for region in module:
            neighbs = [ re.sub('[-\._]', '',x) for x in self.nodes[region] ] 
            region = re.sub('[-\._]','',region)
            print("%s -> {%s};" % (region, ' '.join(neighbs)))
        print("}")

def loadLowMapQRegions(low_mapq_annos):
    low_mapq_regions = {}
    low_mapq_annos = open(low_mapq_annos, 'r')
    for line in low_mapq_annos:
        if line.startswith("#") or line.startswith("chrom"): continue
        toks = line.strip().split('\t')
        region_id = toks[3]
        chrom = toks[6]
        start = int(toks[7])
        end = int(toks[8])
        low_mapq_regions[region_id] = (chrom, start, end)
    low_mapq_annos.close()
    
    return low_mapq_regions

def main(low_mapq_annos, blat_results_file, realign_file, align_to_file, camo_bed):
    regions = loadLowMapQRegions(low_mapq_annos)
    mapping = Graph()
    blat_results = open(blat_results_file, 'r')
    camo_bed = open(camo_bed, 'w')
    for line in blat_results:
        toks = line.strip().split('\t')
        chrom = toks[0]
        start = int(toks[1])
        end = int(toks[2])
        query_name = toks[3]
        target_name = toks[8]

        if target_name != query_name:
            mapping.addEdge(query_name, target_name)
            camo_bed.write("%s\t%d\t%d\n" % (chrom, start, end))
    blat_results.close()

    #modules = mapping.DFS()
    camo_region_realign = open(realign_file, 'w')
    camo_region_alignto = open(align_to_file, 'w')
    nodes = list(mapping.nodes.keys())
    sorted(nodes, key = lambda region_id : len(mapping.nodes[region_id]))
    masked = defaultdict(int)
    for region_id in nodes:
        if masked[region_id] == 0:
            masked[region_id] = 1
            group_ids = [region_id]
            group_pos = [regions[region_id]]
            other_ids = []
            for neighbor_id in mapping.nodes[region_id]:
                if masked[neighbor_id] == 0:
                    group_pos.append(regions[neighbor_id])
                    group_ids.append(neighbor_id)
                    masked[neighbor_id] = 1
                else:
                    other_ids.append(neighbor_id)
            group_ids.extend(other_ids)
            camo_region_alignto.write("%s\t%d\t%d\t" % group_pos[0]) 
            camo_region_alignto.write("%s\t%d\n" % (';'.join(group_ids), len(group_ids)))
            for i in range(len(group_pos)):
                camo_region_realign.write("%s\t%d\t%d\t" % group_pos[i]) 
                camo_region_realign.write("%s\t%d\n" % (';'.join(group_ids[i:] + group_ids[:i]), len(group_ids)))
    camo_region_realign.close()
    camo_region_alignto.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])	
