import pickle
import numpy as np
import ascii_info
import windows_directories
import galcentricutils

# Instantiate Comparitor
compclust = galcentricutils.compclust()

# Define list of clusters to work
saveids = ascii_info.duplimonte_saveids

# Load in the "average" clustering (maindata)
with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    mean_cluster = pickle.load(f)

# Get the current "minimum_excess_index" for unique labelling of the cluster excess
current_excess_value = np.max(mean_cluster) + 1

# Work through all the savids
for saveid in saveids:

    # Load in the clustering
    with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup +
              "\\" + saveid + ".cluster.txt", 'rb') as f:
        current_cluster = pickle.load(f)

    # Compare them
    remap = compclust.compclust_multilabel_julia(mean_cluster, current_cluster, current_excess_value)

    # Set the new excess value
    current_excess_value = np.max(remap) + 1

    # Save the remap
    with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup +
              "\\" + saveid + ".remap-cluster.txt", 'wb') as f:
        current_cluster = pickle.dump(obj=remap, file=f)