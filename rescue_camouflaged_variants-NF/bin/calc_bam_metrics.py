#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

#########################
# Full dark-region file #
#########################
colnames = ['chrom', 'start', 'end', 'MapQBelowThreshold', 'depth', 'percMapQBelowThreshold']
all_dark_regions = pd.read_csv(sys.argv[1], names=colnames, delimiter="\t", dtype={"chrom":str})

# get mean, std, and median depth across the entire genome
mean_depth = all_dark_regions["depth"].mean()
stdev_depth = all_dark_regions["depth"].std()
median_depth = all_dark_regions["depth"].median()



#################################
# Merged low-depth regions file #
#################################

# The mean and median depths in this file are for a given dark region
colnames = ['chrom', 'start', 'end', 'mean depth', 'median depth']
low_depth_regions = pd.read_csv(sys.argv[2], names=colnames, delimiter="\t", dtype={"chrom":str})

# Get number of low-depth regions. This will be total lines in the
# low-depth file, where each line is a dark region. In this case,
# a region is 'dark' if it is <= MIN_DEPTH as defined in 
# combine_DRF_output.py.
n_low_depth_regions = len(low_depth_regions.index)

# Get cumulative size of all dark regions
low_depth_regions['region size'] = low_depth_regions['end'] - low_depth_regions['start']
cum_low_depth_region_size = low_depth_regions['region size'].sum()



################################
# Merged low-mapq regions file #
################################
colnames = ['chrom', 'start', 'end', 'mean depth', 'median depth']
low_mapq_regions = pd.read_csv(sys.argv[3], names=colnames, delimiter="\t", dtype={"chrom":str})

n_low_mapq_regions = len(low_mapq_regions.index)

low_mapq_regions['region size'] = low_mapq_regions['end'] - low_mapq_regions['start']
cum_low_mapq_region_size = low_mapq_regions['region size'].sum()

print("mean_depth\tstdev_depth\tmedian_depth\tn_low_depth_regions\tcum_low_depth_region_size\tn_low_mapq_regions\tcum_low_mapq_region_size")

print("%f\t%f\t%f\t%f\t%f\t%f\t%f" %
        (
            mean_depth,
            stdev_depth,
            median_depth,
            n_low_depth_regions,
            cum_low_depth_region_size,
            n_low_mapq_regions,
            cum_low_mapq_region_size,
    ))
