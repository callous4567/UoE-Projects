import pickle
import numpy as np
import ascii_info_new
import galcentricutils_new
import graphutils_new
import hdfutils
import windows_directories_new

# Instantiate
compclust = galcentricutils_new.compclust()
# Define list of clusters to work
saveids = ascii_info_new.duplimonte_saveids
# Placeholder for all clusterings
clusterings = []
# Work through all the savids
for saveid in saveids:
    # Load the remap
    with open(windows_directories_new.duplimontedir + "\\" + ascii_info_new.fullgroup +
              "\\" + saveid + "_flatfork_.cluster.txt", 'rb') as f:
        current_cluster = pickle.load(file=f)
        clusterings.append(current_cluster)

# Array-up
clustray = np.array(clusterings)

# Work it. probable_clust gives the most probabilistic cluster designation, probability is the probability of it
mapped, memberships = compclust.compclust_multipercentage(clustray, maximum_cluster=np.max(clustray))
# Save it
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.flatfork_asciiname)
writer.write_table(ascii_info_new.fullgroup, "percent_table", mapped)
writer.write_table(ascii_info_new.fullgroup, "total_percent_table", memberships)

# Check the fraction of clusters that actually made it, and save it
viable_fraction = len(clusterings)/len(saveids)
writer.write(ascii_info_new.fullgroup, "viable_cluster_artificial_angular_fraction", viable_fraction)

# Test plot in L-space
probable_clust = mapped['probable_clust']
panda = writer.read_df(ascii_info_new.fullgroup,ascii_info_new.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)
numclust = galcentricutils_new.compclust().nclust_get(probable_clust)
savedexdir = "\\clustered\\flatfork_fullgroup_probable_clust"
data = np.array(data)
graphutils_new.threed_graph().kmeans_L_array(data, probable_clust, savedexdir, browser=False,  outliers=False)

# Spatially, too
posdata = np.array([panda['x'],panda['y'],panda['z']]).T
graphutils_new.threed_graph().xmeans_L_array(posdata, probable_clust, savedexdir + "_spatial", browser=False,  outliers=False)

# Test out the vasiliev graphing for this set...
table = writer.read_table(ascii_info_new.fullgroup,ascii_info_new.fullset)
for clust_id in range(0, numclust):
    graphutils_new.spec_graph().clust_radec(table,
                                        probable_clust,
                                        cluster_id=clust_id,
                                        savedexdir="flatfork_memberpercent_" + str(clust_id) + "_test",
                                        vasiliev=False,
                                        lb=True)