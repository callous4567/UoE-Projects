import pickle
import numpy as np
from matplotlib import pyplot as plt
from munkres import Munkres
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.cluster import contingency_matrix
import ascii_info
import seaborn as sns
import time
import galcentricutils
import hdfutils
import graphutils
import windows_directories


panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                             ascii_info.set_raw)

# Load in the lists as a list of lists! Clustering arrays.
arrays = [windows_directories.duplimontedir + "\\" + ascii_info.fullgroup + "\\" + d + ".cluster.txt" \
          for d in ascii_info.duplimonte_saveids]
nclusts = []
for num, array in enumerate(arrays):
    with open(array, 'rb') as f:
        arrays[num] = pickle.load(f)
        owo = galcentricutils.compclust().nclust_get_complete(arrays[num])
        nclusts.append(galcentricutils.compclust().nclust_get(arrays[num]))

# Array Choices
numone, numtwo = 0, 1

new_array = galcentricutils.compclust().compclust_multi(arrays[numone], arrays[numtwo])
#for one,two in zip(arrays[numone], new_array):
    #print(one, two)


# Load in the relevant data arrays, also.
L_arrays = [windows_directories.duplimontedir + "\\" + ascii_info.fullgroup + "\\" + d + ".txt" \
          for d in ascii_info.duplimonte_saveids]
for num, array in enumerate(L_arrays):
    with open(array, 'rb') as f:
        L_arrays[num] = pickle.load(f)
print(L_arrays[0], type(L_arrays[0]))

graphutils.threed_graph().kmeans_L_array_pair([L_arrays[numone],L_arrays[numtwo]], [arrays[numone],new_array], False, False, False)