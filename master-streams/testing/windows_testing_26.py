import pickle
import numpy as np
from pandas import DataFrame

import ascii_info
import hdfutils
import windows_directories

# Maindata clustering
with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(file=f)

# Generate a table of memberships
prelimembers = [0 for d in np.arange(-1, np.max(clustered) + 1)]
for cluster in np.arange(-1, np.max(clustered) + 1, 1):
    for star in clustered:
        if star == cluster:
            prelimembers[cluster + 1] += 1
prelimembers = np.array(prelimembers)

# Obtain the monte carlo realisation membership
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
membership_monte = writer.read_table(ascii_info.fullgroup, "percent_table")
percenttable_members = [0 for d in np.arange(-1, np.max(clustered) + 1)]
for cluster in np.arange(-1, np.max(clustered) + 1, 1):
    for star in membership_monte['probable_clust']:
        if star == cluster:
            percenttable_members[cluster + 1] += 1
percenttable_members = np.array(percenttable_members)

# And the membership with greatcircles
membership_table_monte_gc = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")
greattable_members = [0 for d in np.arange(-1, np.max(clustered) + 1)]
for cluster in np.arange(-1, np.max(clustered) + 1, 1):
    for star in membership_table_monte_gc['greatcircle_probable_clust']:
        if star == cluster:
            greattable_members[cluster + 1] += 1
greattable_members = np.array(greattable_members)

# Set data array up
data_array = np.zeros((np.max(clustered) + 2, 4))
data_array[:,0] = np.arange(-1, np.max(clustered) + 1)
data_array[:,1] = prelimembers
data_array[:,2] = percenttable_members
data_array[:,3] = greattable_members

dataframe = DataFrame(data_array, columns=["Cluster", "Preliminary", "Monte", "Monte GCCs"])
dataframe.to_latex(caption="Cluster memberships determined for the dataset: for the preliminary clustering, quasi-soft clustering, and quasi-soft refined clustering. The cluster index -1 corresponds to the stars labelled as ``noise'' by \verb|hdbscan|.",
                   buf=windows_directories.imgdir + "\\memberships.tex",
                   label="tab:memberships",
                   float_format="%.0d",
                   index=False)
