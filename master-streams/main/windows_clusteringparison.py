import numpy as np
import galcentricutils
import ascii_info

# Visualize number of clusters caught for each Monte-Carloed angular momentum set (histogram included as side-bar.)
galcentricutils.compclust().graph_nclust(ascii_info.bhb)
galcentricutils.compclust().graph_nclust(ascii_info.kgiants)
galcentricutils.compclust().graph_nclust(ascii_info.gcs)
galcentricutils.compclust().graph_nclust(ascii_info.lamostk)
galcentricutils.compclust().graph_nclust(ascii_info.fullgroup)


"""
# Compare two clusterings- test data
import graphutils
import windows_directories
import pickle

with open(windows_directories.duplimontedir + "\\full_raw\\L_0.cluster.txt", 'rb') as f:
    clust1 = pickle.load(f)

with open(windows_directories.duplimontedir + "\\full_raw\\L_100.cluster.txt", 'rb') as f:
    clust2 = pickle.load(f)

with open(windows_directories.duplimontedir + "\\full_raw\\L_0.txt", 'rb') as f:
    data1 = pickle.load(f)

with open(windows_directories.duplimontedir + "\\full_raw\\L_100.txt", 'rb') as f:
    data2 = pickle.load(f)

clust2_relabel = galcentricutils.compclust().compclust(clust1,clust2)
graphutils.threed_graph().kmeans_L_array_pair([data1,data2],[clust1,clust2_relabel],False,True,outliers=False) """