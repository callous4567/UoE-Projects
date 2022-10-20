import graphutils_new
import pickle
import numpy as np
import ascii_info_new
import galcentricutils_new
import hdfutils
import windows_directories_new

# Preliminary Clustering for Greatcount
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.fullset)
clustered = np.array(table['prelim_clust'], int)
numclust = galcentricutils_new.compclust().nclust_get(clustered)

# Get the unique clusters in this
clusters_to_try = []
for clust in clustered:
    if clust not in clusters_to_try:
        clusters_to_try.append(clust)

# Define the "threshold" for a cluster to be considered (the number of populated clusts: see cluster_comparison.)
clust_threshold = 4
# generally speaking, this seems suitable- each clustering, even with a greatcircle, has half
# the usual number available.

# Instantiate
compclust = galcentricutils_new.compclust()

# Grab saveids
saveids = ascii_info_new.duplimonte_LXYZ_saveids

# Placeholder for all clusterings (stored as tables.)
clusterings = []

run = True
if run == True:

    # For each cluster, create a percenttable.
    for cluster in clusters_to_try:

        # List of all clusterings
        clustering_list = []

        # Load in the clusterings from all saveids for this cluster
        for saveid in saveids:

            # Load the remap
            with open(windows_directories_new.duplimontedir + "\\" + ascii_info_new.fullgroup +
                      "\\" + saveid + (".remap-cluster_{0:.0f}.txt").format(cluster), 'rb') as f:
                current_cluster = pickle.load(file=f)

                # Verify the number of actually populated clusters in current_cluster
                nclust = compclust.nclust_get_complete(current_cluster)

                # Only memberpercent it if it satisfies n > nclust
                if nclust[0] > clust_threshold:
                    clustering_list.append(current_cluster)

        # Array it up
        clustering_list = np.array(clustering_list)

        # Generate memberships/etc. First one has all percentages, second one just has individual ones.
        mapped, memberships = compclust.compclust_multipercentage(clustering_list,
                                                                  maximum_cluster=np.max(clustering_list))

        # Append to the table for this particular cluster
        clusterings.append(mapped)

        # Save the table, too.
        writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
        writer.write_table(ascii_info_new.fullgroup, ("percent_table_{0:.0f}").format(cluster), mapped)
        writer.write_table(ascii_info_new.fullgroup, ("total_percent_table_{0:.0f}").format(cluster), memberships)

else:
    # Load in the data
    for cluster in clusters_to_try:
        # Read the table
        writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
        mapped = writer.read_table(ascii_info_new.fullgroup, ("percent_table_{0:.0f}").format(cluster))

        # Append
        clusterings.append(mapped)

# Load in the "mean" percenttable and map
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
meanmap = writer.read_table(ascii_info_new.fullgroup, "percent_table")

# Get the final membership table/map.
meanmap, meanmemberships = compclust.percenttable_stack(meanmap, clusterings, clusters_to_try)
writer.write_table(ascii_info_new.fullgroup, "percent_table_greatfitted", meanmap), \
writer.write_table(ascii_info_new.fullgroup, "total_percent_table_greatfitted", meanmemberships)

# Create final plots under "clustered"
numclust = galcentricutils_new.compclust().nclust_get(meanmap['probable_clust'])
savedexdir = "\\clustered\\fullgroup_greatfit_clustered_" + str(numclust) + "_clusters_found"
panda = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_df(ascii_info_new.fullgroup,
                                                           ascii_info_new.panda_raw)
data = np.array(panda['vec_L'])
data = list(data)
graphutils_new.threed_graph().kmeans_L_array(data, meanmap['probable_clust'], savedexdir, browser=False, outliers=False)

# Test out the vasiliev graphing for this set...
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.fullset)
for clust_id in range(0, numclust):
    graphutils_new.spec_graph().clust_radec(table, meanmap['probable_clust'], cluster_id=clust_id,
                                        savedexdir="greatfit" + "_" + str(clust_id) + "_clusttest")
