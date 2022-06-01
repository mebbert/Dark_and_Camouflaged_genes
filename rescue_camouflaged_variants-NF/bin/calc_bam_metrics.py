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
all_mean_depth = all_dark_regions["depth"].mean()
all_stdev_depth = all_dark_regions["depth"].std()
all_median_depth = all_dark_regions["depth"].median()



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


low_depth_mean_depth = low_depth_regions["mean depth"].mean()
low_depth_stdev_depth = low_depth_regions["mean depth"].std()
low_depth_median_depth = low_depth_regions["mean depth"].median()


################################
# Merged low-mapq regions file #
################################
colnames = ['chrom', 'start', 'end', 'mean depth', 'median depth']
low_mapq_regions = pd.read_csv(sys.argv[3], names=colnames, delimiter="\t", dtype={"chrom":str})

n_low_mapq_regions = len(low_mapq_regions.index)

low_mapq_regions['region size'] = low_mapq_regions['end'] - low_mapq_regions['start']
cum_low_mapq_region_size = low_mapq_regions['region size'].sum()

low_mapq_mean_depth = low_mapq_regions["mean depth"].mean()
low_mapq_stdev_depth = low_mapq_regions["mean depth"].std()
low_mapq_median_depth = low_mapq_regions["mean depth"].median()



###############################
# Write stats to file         #
###############################

print("all_mean_depth\tall_stdev_depth\tall_median_depth\tn_low_depth_regions\tcum_low_depth_region_size\tlow_depth_mean_depth\tlow_depth_stdev_depth\tlow_depth_median_depth\tn_low_mapq_regions\tcum_low_mapq_region_size\tlow_mapq_mean_depth\tlow_mapq_stdev_depth\tlow_mapq_median_depth")

print("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" %
        (
            all_mean_depth,
            all_stdev_depth,
            all_median_depth,
            n_low_depth_regions,
            cum_low_depth_region_size,
            low_depth_mean_depth,
            low_depth_stdev_depth,
            low_depth_median_depth,
            n_low_mapq_regions,
            cum_low_mapq_region_size,
            low_mapq_mean_depth,
            low_mapq_stdev_depth,
            low_mapq_median_depth,
    ))
