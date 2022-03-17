import pickle
import numpy as np
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories

"""
Candidate function to handle great circle and clustering.
"""
# Clustering parameters/etc
arrayinfominpars = []
group = ascii_info.fullgroup
minpar = ascii_info.fulldata_minpar

# Set up data and run the clustering
panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)
clustered = galcentricutils.cluster3d().listhdbs(data, minpar)

# Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
numclust = galcentricutils.compclust().nclust_get(clustered)
savedexdir = "\\clustered\\fullgroup_clustered_" + str(numclust) + "_clusters_found"
data = np.array(data)
graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=True,  outliers=False)


# Test out the vasiliev graphing for this set...
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
for clust_id in range(0, numclust):
    graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id, savedexdir=str(clust_id) + "_clusttest")