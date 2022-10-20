import pickle
import numpy as np
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories

# Instantiate
compclust = galcentricutils.compclust()
# Define list of clusters to work
saveids = ascii_info.duplimonte_saveids
# Placeholder for all clusterings
clusterings = []
# Work through all the savids
for saveid in saveids:
    # Load the remap
    with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup +
              "\\" + saveid + "_finetune.remap-cluster.txt", 'rb') as f:
        current_cluster = pickle.load(file=f)
        clusterings.append(current_cluster)

# Array-up
clusterings = np.array(clusterings)

# Work it. probable_clust gives the most probabilistic cluster designation, probability is the probability of it
mapped, memberships = compclust.compclust_multipercentage(clusterings, maximum_cluster=np.max(clusterings))

# Save it
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.finetune_asciiname)
writer.write_table(ascii_info.fullgroup, "percent_table", mapped)
writer.write_table(ascii_info.fullgroup, "total_percent_table", memberships)

# Check the fraction of clusters that actually made it, and save it
viable_fraction = len(clusterings)/len(saveids)
writer.write(ascii_info.fullgroup, "viable_cluster_artificial_angular_fraction", viable_fraction)

# Test plot in L-space
probable_clust = mapped['probable_clust']
panda = writer.read_df(ascii_info.fullgroup,ascii_info.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)
numclust = galcentricutils.compclust().nclust_get(probable_clust)
savedexdir = "\\clustered\\finetune_fullgroup_probable_clust"
data = np.array(data)
graphutils.threed_graph().kmeans_L_array(data, probable_clust, savedexdir, browser=False,  outliers=False)

# Spatially, too
posdata = np.array([panda['x'],panda['y'],panda['z']]).T
graphutils.threed_graph().xmeans_L_array(posdata, probable_clust, savedexdir + "_spatial", browser=False,  outliers=False)

# Test out the vasiliev graphing for this set...
table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)
for clust_id in range(0, numclust):
    graphutils.spec_graph().clust_radec(table,
                                        probable_clust,
                                        cluster_id=clust_id,
                                        savedexdir="finetune_memberpercent_" + str(clust_id) + "_test",
                                        vasiliev=False,
                                        flatfork="finetune")