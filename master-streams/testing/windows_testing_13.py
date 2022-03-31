import itertools
import pickle

import ascii_info
import galcentricutils
import hdfutils
import windows_directories
from energistics import orbigistics
from galpy.util import bovy_coords
import numpy as np

from graphutils import twod_graph

# Load in the "mean" percenttable and map
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
membership_table = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")
table = writer.read_table(ascii_info.fullgroup, ascii_info.fullset)
table = table[[True if d == 8 else False for d in membership_table['probable_clust']]]
table = galcentricutils.angular().get_polar(table)
print(table['r'])