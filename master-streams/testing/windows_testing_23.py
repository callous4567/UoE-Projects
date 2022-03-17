import copy
import itertools
import numpy as np

# Load in the "mean" percenttable and map
import ascii_info
import hdfutils
import windows_directories
from galcentricutils import compclust

writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
meanmap = writer.read_table(ascii_info.fullgroup, "percent_table")
for i in range(len(meanmap)):
    if meanmap['probable_clust'][i] == 1:
        print(i)