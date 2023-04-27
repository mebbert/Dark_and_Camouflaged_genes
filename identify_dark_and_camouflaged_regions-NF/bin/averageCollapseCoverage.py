#!/usr/bin/env python3
import sys
import numpy as np



for line in sys.stdin:
    toks = line.strip().split('\t')
    pos = toks[0:3]
    collapsedTotalDepth=[float(x) for x in toks[3].split(",")]
    meanDepth = np.mean(collapsedTotalDepth)
    medianDepth = np.median(collapsedTotalDepth)

    line="\t".join(pos)+"\t"+str(meanDepth)+"\t"+str(medianDepth)+"\n"
    sys.stdout.write(line)
