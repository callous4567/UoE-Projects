import multiprocessing
import warnings
import ascii_info, windows_directories
import numpy as np
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

# Preliminary Clustering for Greatcount
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
clustered = np.array(table['prelim_clust'], int)
numclust = galcentricutils.compclust().nclust_get(clustered)

# Get the unique clusters in this
clusters_to_try = []
for clust in clustered:
    if clust not in clusters_to_try:
        clusters_to_try.append(clust)

gcc_widths = [30 for d in clusters_to_try] # roughly twice the width of the stream, I would say.
arrayinfominpars = []

for saveid in ascii_info.duplimonte_LXYZ_saveids:
    arrayinfominpars.append([[group,saveid],minpar,clusters_to_try, gcc_widths])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.do_hdbscan_greatfit, arrayinfominpars)
    pool.close()


