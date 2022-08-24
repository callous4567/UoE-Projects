import itertools
import pickle

import ascii_info_new
import galcentricutils_new
import hdfutils
import windows_directories_new
from energistics_new import orbigistics
from galpy.util import bovy_coords
import numpy as np

from graphutils_new import twod_graph

# Load in the "mean" percenttable and map
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
membership_table = writer.read_table(ascii_info_new.fullgroup, "percent_table_greatfitted")
table = writer.read_table(ascii_info_new.fullgroup, ascii_info_new.fullset)
table = table[[True if d == 8 else False for d in membership_table['probable_clust']]]
table = galcentricutils_new.angular().get_polar(table)
print(table['r'])