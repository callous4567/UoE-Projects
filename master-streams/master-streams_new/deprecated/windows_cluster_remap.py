#windows_directories_new.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".txt", 'rb'
#group = ascii_info_new.fullgroup
#minpar = ascii_info_new.fulldata_minpar
import pickle
import graphutils_new
import numpy as np
import galcentricutils_new
import ascii_info_new
import hdfutils
import windows_directories_new

# Load the "main" clustering of the mean data
with open(windows_directories_new.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(f)


# Set up data and run the clustering
panda = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_df(ascii_info_new.fullgroup,
                                                           ascii_info_new.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)

# Set up a list of all objects to "match" against
#clusterings = [d + ".cluster.txt" for d in ascii_info_new.duplimonte_saveids]
#
# List to iterate
#ref_and_clust = []
#for clustering in clusterings:
#    ref_and_clust.append([clustered, clustering])
#
# Test
with open(windows_directories_new.duplimontedir + "\\" + ascii_info_new.fullgroup + "\\" + "L_4.cluster.txt", 'rb') as f:
    cluster_test = pickle.load(f)
remap_test = galcentricutils_new.compclust().compclust_multilabel(cluster_test, clustered, np.max(clustered)) # clustered,

# Grab test data
with open(windows_directories_new.duplimontedir + "\\" + ascii_info_new.fullgroup + "\\" + "L_4.txt", 'rb') as f:
    remap_data = pickle.load(f)

graphutils_new.threed_graph().kmeans_L_array_pair([data, remap_data], [clustered, remap_test], False, True, False)