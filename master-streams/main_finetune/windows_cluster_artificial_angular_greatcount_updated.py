import multiprocessing
import warnings
import ascii_info, windows_directories
import numpy as np
import hdfutils
import galcentricutils
import windows_multiprocessing

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

writer = hdfutils.hdf5_writer(windows_directories.datadir,
                              ascii_info.finetune_asciiname)
if __name__ == "__main__":

    # Set up variables for clusterings
    group = ascii_info.fullgroup
    min_clust_size = ascii_info.fulldata_minpar[0] # min_samples from greattable

    # Grab the finetune table *and* greattable
    table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)

    greattable = writer.read_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles"
    )

    # Add a wince to the distance max_dist to ensure boundaries aren't occluded
    greattable['max_dist'] += 1

    clustered = np.array(table['prelim_clust'], int)
    numclust = galcentricutils.compclust().nclust_get(clustered)

    # Choose clusts to try greatcount clustering for
    clusters_to_try = []
    for clust in list(set(clustered)):

        if clust not in clusters_to_try:

            if clust != -1:

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if len(np.where(clustered==clust)[0])/len(clustered) < 0.1:

                    clusters_to_try.append(clust)
    greattable_indices = [
        np.where(greattable['clust_to_try']==clust)[0][0] for clust in clusters_to_try
    ]
    best_minpars = greattable['best_samples'][greattable_indices]
    best_minpars = [[ascii_info.fulldata_minpar[0], min_samples] for min_samples in best_minpars]
    widths = greattable['max_dist'][greattable_indices]
    thetas = greattable['theta'][greattable_indices]
    phis = greattable['phi'][greattable_indices]
    greatpars = np.array([widths,thetas,phis]).T
    arrayinfominpars = []
    for saveid in ascii_info.duplimonte_LXYZ_saveids:
        arrayinfominpars.append([[group,saveid],best_minpars,clusters_to_try, greatpars])

    # Regularly map/pool :)
    pool = multiprocessing.Pool(4)
    results = pool.map(windows_multiprocessing.do_updated_hdbscan_greatfit_finetune, arrayinfominpars)
    pool.close()


