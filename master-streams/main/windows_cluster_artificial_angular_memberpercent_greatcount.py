import graphutils
from ascii_info import clusters_to_cluster
import pickle
import numpy as np
import ascii_info
import galcentricutils
import hdfutils
import windows_directories

# Instantiate
compclust = galcentricutils.compclust()

# Grab saveids
saveids = ascii_info.duplimonte_LXYZ_saveids

# Placeholder for all clusterings (stored as tables.)
clusterings = []

run = True
if run == True:
    # For each cluster, create a percenttable.
    for cluster in clusters_to_cluster:

        # List of all clusterings
        clustering_list = []

        # Load in the clusterings from all saveids for this cluster
        for saveid in saveids:

            # Load the remap
            with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup +
                      "\\" + saveid + (".remap-cluster_{0:.0f}.txt").format(cluster), 'rb') as f:
                clustering_list.append(pickle.load(file=f))

        # Array it up
        clustering_list = np.array(clustering_list)

        # Generate memberships/etc. First one has all percentages, second one just has individual ones.
        mapped, memberships = compclust.compclust_multipercentage(clustering_list, maximum_cluster=np.max(clustering_list))

        # Append to the table for this particular cluster
        clusterings.append(mapped)

        # Save the table, too.
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        writer.write_table(ascii_info.fullgroup, ("percent_table_{0:.0f}").format(cluster), mapped)
        writer.write_table(ascii_info.fullgroup, ("total_percent_table_{0:.0f}").format(cluster), memberships)

else:
    # Load in the data
    for cluster in clusters_to_cluster:

        # Read the table
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        mapped = writer.read_table(ascii_info.fullgroup, ("percent_table_{0:.0f}").format(cluster))

        # Append
        clusterings.append(mapped)

# Load in the "mean" percenttable and map
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
meanmap = writer.read_table(ascii_info.fullgroup, "percent_table")

# Get the final membership table/map.
meanmap, meanmemberships = compclust.percenttable_stack(meanmap, clusterings, clusters_to_cluster)
writer.write_table(ascii_info.fullgroup, "percent_table_greatfitted", meanmap), \
writer.write_table(ascii_info.fullgroup, "total_percent_table_greatfitted", meanmemberships)

# Create final plots under "clustered"
numclust = galcentricutils.compclust().nclust_get(meanmap['probable_clust'])
savedexdir = "\\clustered\\fullgroup_greatfit_clustered_" + str(numclust) + "_clusters_found"
panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)
graphutils.threed_graph().kmeans_L_array(data, meanmap['probable_clust'], savedexdir, browser=False,  outliers=False)


# Test out the vasiliev graphing for this set...
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
for clust_id in range(0, numclust):
    graphutils.spec_graph().clust_radec(table, meanmap['probable_clust'], cluster_id=clust_id, savedexdir="greatfit" + "_" + str(clust_id) + "_clusttest")

