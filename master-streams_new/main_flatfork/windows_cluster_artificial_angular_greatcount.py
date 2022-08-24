import multiprocessing
import warnings
import ascii_info_new, windows_directories_new
import pickle
import numpy as np
import graphutils_new
import hdfutils
import galcentricutils_new
import windows_multiprocessing_new

# TODO: Note that this is the old version without fine-tuning of parameters.

"""
Second round of clustering: cluster, for each greatcircle, all the monte-carlo-produced data. 
- Pass the combination of minpar, cluster_to_cluster, and saveid
- Load in that saveid monte-carlo in LXYZ 
- For each cluster_to_cluster, do GCC, retaining the length of the table for comparative purposes, then cluster it
- Compare the clustering to the mean clustering (thus matching to the original cluster_to_cluster) and relabel
- Save this relabelling
"""

# Suppress Loky Warnings (don't ask don't tell- damned annoying.)
warnings.filterwarnings("ignore",
                        message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

if __name__ == "__main__":

    # Set up variables for clusterings
    group = ascii_info_new.fullgroup
    minpar = ascii_info_new.fulldata_minpar # for hdbscan. use the oldminpar for these sets, I reckon.

    # Grab the flatfork table *and* greattable
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                           ascii_info_new.fullset)

    greattable = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                      ascii_info_new.flatfork_asciiname).read_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles")

    # Add a wince to the distance max_dist to ensure boundaries aren't occluded
    greattable['max_dist'] += 1

    clustered = np.array(table['prelim_clust'], int)
    numclust = galcentricutils_new.compclust().nclust_get(clustered)

    # Get the unique clusters in this
    clusters_to_try = []
    for clust in clustered:
        if clust not in clusters_to_try:
            clusters_to_try.append(clust)

    gcc_widths = [30 for d in clusters_to_try] # roughly twice the width of the stream, I would say.
    arrayinfominpars = []

    for saveid in ascii_info_new.duplimonte_LXYZ_saveids:
        arrayinfominpars.append([[group,saveid],minpar,clusters_to_try, gcc_widths])

    # Regularly map/pool :)
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing_new.do_hdbscan_greatfit, arrayinfominpars)
    pool.close()


