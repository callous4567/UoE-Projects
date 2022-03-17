import multiprocessing
import warnings
import ascii_info, windows_directories
import pickle
import numpy as np
import graphutils
import hdfutils
import galcentricutils
import windows_multiprocessing

"""
Second round of clustering: cluster, for each greatcircle, all the monte-carlo-produced data. 
- Pass the combination of minpar, cluster_to_cluster, and saveid
- Load in that saveid monte-carlo in LXYZ 
- For each cluster_to_cluster, do GCC, retaining the length of the table for comparative purposes, then cluster it
- Compare the clustering to the mean clustering (thus matching to the original cluster_to_cluster) and relabel
- Save this relabelling
"""

# Suppress Loky Warnings
warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

# Set up variables for clusterings
group = ascii_info.fullgroup
minpar = ascii_info.fulldata_minpar # for hdbscan. use the oldminpar for these sets, I reckon.
clusters_to_cluster = ascii_info.clusters_to_cluster
gcc_widths = ascii_info.gcc_widths # roughly twice the width of the stream, I would say.
arrayinfominpars = []

for saveid in ascii_info.duplimonte_LXYZ_saveids:
    arrayinfominpars.append([[group,saveid],minpar,clusters_to_cluster, gcc_widths])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.do_hdbscan_greatfit, arrayinfominpars)
    pool.close()


