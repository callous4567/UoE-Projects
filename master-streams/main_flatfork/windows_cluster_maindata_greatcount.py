import pickle
import numpy as np
from astropy.table import Table
from matplotlib import tri
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
from energistics import energistics, energistics_manual, orbigistics
from galcentricutils import greatcount
import matplotlib.pyplot as plt


just_plot = False

if just_plot == False and __name__ == "__main__":

    # Preliminary Clustering for Greatcount
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.flatfork_asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)
    table = galcentricutils.angular().get_polar(table)
    clustered = np.array(table['prelim_clust'], int)
    numclust = galcentricutils.compclust().nclust_get(clustered)

    # Get the unique clusters in this (save noise)
    clusters_to_try = []
    for clust in list(set(clustered)):

        if clust not in clusters_to_try:

            if clust != -1:

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if len(np.where(clustered==clust)[0])/len(clustered) < 0.1:

                    clusters_to_try.append(clust)

    # Generate the grid of great circles to try for the "first pass."
    grid_edge_length = 800 # number of GCC poles to generate
    resolution = 800 # number of datapoints to generate in GCC
    pole_thetas, pole_phis, pole_indices = greatcount().fullgridgen(np.rint(grid_edge_length**2, casting='unsafe'))
    pole_thetas, pole_phis = np.degrees(pole_thetas), np.degrees(pole_phis)

    # Generate run class
    greatfit = galcentricutils.greatfit()

    # Do this method if you're dealing with "real_dist=False" and just least-squaring.
    """
    real_dist will specify the distance for which to maximize enclosed points to the great circle: 
    see the implementation 
    In degrees!
    """
    real_dist = 5

    # Make the zip required for us
    zipped = []
    for which_try in clusters_to_try:

        table_included = table[[True if d == which_try else False for d in clustered]]
        clustered_included = table_included['prelim_clust']
        zip = [real_dist, which_try, table_included, clustered_included, pole_thetas, pole_phis, resolution]
        zipped.append(zip)

    # Map
    import multiprocessing, windows_multiprocessing
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.flatfork_greatcircle_optimize, zipped)
    pool.close()

    # Format results from results
    poles_found, poles_stdevs, poles_maxdists = [],[],[]
    for result in results:
        poles_found.append(result[0])
        poles_stdevs.append(result[1])
        poles_maxdists.append(result[2])

    # Save the best_circles to a table (non-memberpercent).
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.flatfork_asciiname)
    table = Table()
    table['clust_to_try'] = clusters_to_try
    thetapoles, phipoles = np.array(poles_found).T
    table['theta'] = thetapoles
    table['phi'] = phipoles
    table['stdev_dist'] = poles_stdevs
    table['max_dist'] = poles_maxdists
    writer.write_table("greatcircles_nonmemberpercent_preliminary", "greatcircles", table)

if just_plot == True and __name__ == "__main__":

    # Preliminary Clustering for Greatcount
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.flatfork_asciiname).read_table(ascii_info.fullgroup,
                                                                           ascii_info.fullset)
    table = galcentricutils.angular().get_polar(table)
    clustered = np.array(table['prelim_clust'], int)
    numclust = galcentricutils.compclust().nclust_get(clustered)
    greatcircle_table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.flatfork_asciiname).read_table("greatcircles_nonmemberpercent_preliminary",
                                                                           "greatcircles")

    # Get the unique clusters in this (save noise)
    clusters_to_try = []
    for clust in list(set(clustered)):

        if clust not in clusters_to_try:

            if clust != -1:

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if len(np.where(clustered==clust)[0])/len(clustered) < 0.1:

                    clusters_to_try.append(clust)

    # Generate the grid of great circles to try for the "first pass."
    grid_edge_length = 400  # number of GCC poles to generate
    resolution = 400  # number of datapoints to generate in GCC

    # Generate run class
    greatfit = galcentricutils.greatfit()

    # For each individual clustering we have, try to produce plots
    for which_try in clusters_to_try:

        # Demarcate which ones to use
        included = [True if d == which_try else False for d in clustered]
        table_included = table[included]

        # Get the row of the table
        arg = np.where(greatcircle_table['clust_to_try']==which_try)[0][0]
        best_circle = greatcircle_table[arg][['theta','phi']]

        # Go forth and plot the greatcircle against the data, in ra/dec.
        try:
            graphutils.spec_graph().clust_thetaphi(table=table_included, clustering=clustered[included],
                                                   cluster_id=which_try,
                                                   vasiliev=False, savedexdir="flatfork_greatcount_beforeMC_" + str(
                    which_try) + "_thetaphi_greatcircle",
                                                   gcc_thetaphis=greatfit.gcc_gen(1000, *best_circle),
                                                   flatfork=True)
            import matplotlib.pyplot as plt

            plt.close()
            plt.clf()
        except Exception as e:
            print(e)
            pass

