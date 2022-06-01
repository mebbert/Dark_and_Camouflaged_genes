#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import sys

#################################
# Merged low-depth regions file #
#################################

colnames = ["extractionChrom", "extractionStart", "extractionEnd", "extractionAnnotation", "numAnnotations", "chrom", "start", "end", "nMAPQThreshold", "depth", "percentMapqThreshold"]
extraction_intersect = pd.read_csv(sys.argv[1], names=colnames, delimiter="\t")

extraction_intersect["region"] = extraction_intersect["extractionChrom"] + ":" + extraction_intersect["extractionStart"].astype(str) + "-" + extraction_intersect["extractionEnd"].astype(str)

uniqueRegions = extraction_intersect["region"].unique()

output = pd.DataFrame(columns = ["chrom", "start", "end", "region", "annoation", "meanDepth", "medianDepth", "stdevDepth"])

for region in uniqueRegions:
    r = re.split(":|-", region)
    subset = extraction_intersect.where((extraction_intersect["region"] == region) & (extraction_intersect["chrom"] != "."))
    subset.replace('.',np.NaN)
    if len(subset.index) > 0:
        mean=subset["depth"].astype(float).mean()
        median=subset["depth"].astype(float).median()
        stdev=subset["depth"].astype(float).std()

        output.append({"chrom": r[0], "start": r[1], "end": r[2],"regions": region, "annotation": subset["extractionAnnotation"][0], "meanDepth": mean, "medianDepth": median, "stdevDepth": stdev}, ignore_index=True)




output.to_csv(sys.argv[2], sep="\t")

