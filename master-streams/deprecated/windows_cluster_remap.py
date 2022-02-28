#windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".txt", 'rb'
#group = ascii_info.fullgroup
#minpar = ascii_info.fulldata_minpar
import pickle
import graphutils
import numpy as np
import galcentricutils
import ascii_info
import hdfutils
import windows_directories

# Load the "main" clustering of the mean data
with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(f)


# Set up data and run the clustering
panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)

# Set up a list of all objects to "match" against
#clusterings = [d + ".cluster.txt" for d in ascii_info.duplimonte_saveids]
#
# List to iterate
#ref_and_clust = []
#for clustering in clusterings:
#    ref_and_clust.append([clustered, clustering])
#
# Test
with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup + "\\" + "L_4.cluster.txt", 'rb') as f:
    cluster_test = pickle.load(f)
remap_test = galcentricutils.compclust().compclust_multilabel(cluster_test, clustered, np.max(clustered)) # clustered,

# Grab test data
with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup + "\\" + "L_4.txt", 'rb') as f:
    remap_data = pickle.load(f)

graphutils.threed_graph().kmeans_L_array_pair([data, remap_data], [clustered, remap_test], False, True, False)