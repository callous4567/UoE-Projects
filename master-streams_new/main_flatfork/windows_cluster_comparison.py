import numpy as np
import galcentricutils_new
import ascii_info_new

# Visualize number of clusters caught for each Monte-Carloed angular momentum set (histogram included as side-bar.)
#galcentricutils_new.compclust().graph_nclust(ascii_info_new.fullgroup)

# Visualize number of clusters caught for each Monte-Carloed angular momentum set (histogram included as side-bar.)
# Set energy to True to visualize from [Lx Ly Lz E] clustering and Sofie for L LZ E CIRC
galcentricutils_new.compclust().flatfork_graph_nclust(ascii_info_new.fullgroup, energy=False, sofie=False)

# Example illustrating use of pair-comparison tool of graphutils_new. TODO: Find a place in actual documentation for this.
"""
# Compare two clusterings- test data 
import graphutils_new
import windows_directories_new
import pickle

with open(windows_directories_new.duplimontedir + "\\full_raw\\L_0.cluster.txt", 'rb') as f:
    clust1 = pickle.load(f)
    print(np.max(clust1))

with open(windows_directories_new.duplimontedir + "\\full_raw\\L_100.cluster.txt", 'rb') as f:
    clust2 = pickle.load(f)
    print(np.max(clust2))

with open(windows_directories_new.duplimontedir + "\\full_raw\\L_0.txt", 'rb') as f:
    data1 = pickle.load(f)

with open(windows_directories_new.duplimontedir + "\\full_raw\\L_100.txt", 'rb') as f:
    data2 = pickle.load(f)

clust2_relabel = galcentricutils_new.compclust().\
    compclust_multilabel(clust2,clust1,
                         minimum_excess_index=np.max(np.array([np.max(clust1), np.max(clust2)])) + 1)
print(np.max(clust2_relabel))
graphutils_new.threed_graph().kmeans_L_array_pair([data1,data2],[clust1,clust2_relabel],False,True,outliers=False) """