import pandas as pd
import numpy as np
import sys

bases = pd.read_csv(sys.argv[1], delimiter="\t", dtype={"chrom":str})
mean = bases["depth"].mean()
median = bases["depth"].median()
print("%f\t%f" % (mean, median))
