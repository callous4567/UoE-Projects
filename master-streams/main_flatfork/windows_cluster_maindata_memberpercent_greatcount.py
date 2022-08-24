import numpy as np
from astropy.table import Table
import ascii_info
import galcentricutils
import hdfutils
import windows_directories
from galcentricutils import greatcount

# Load the post-initial-montecarlo-clustering (via memberpercent.py.)
"""
This has the benefit of having weeded-out most of the "non-clusterings" via monte-carlo-
most of them will have vanished by this point and been absorbed into other clusters, if they're not the "significant one."
"""

if __name__ == "__main__":

    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    percent_table = writer.read_table(ascii_info.fullgroup, "percent_table")
    clustered = percent_table['probable_clust']

    #with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    #    clustered = pickle.load(file=f)

    # Get the unique clusters in this (save noise)
    clusters_to_try = []
    for clust in list(set(clustered)):

        if clust not in clusters_to_try:

            if clust != -1:

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if len(np.where(clustered==clust)[0])/len(clustered) < 0.1:

                    clusters_to_try.append(clust)


    # Set up data
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)
    table = galcentricutils.angular().get_polar(table)

    # Generate the grid of great circles to try for the "first pass"
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
    results = pool.map(windows_multiprocessing.flatfork_greatcircle_optimize_memberpercent, zipped)
    pool.close()

    # Format results from results
    poles_found, poles_stdevs, poles_maxdists = [],[],[]
    for result in results:
        poles_found.append(result[0])
        poles_stdevs.append(result[1])
        poles_maxdists.append(result[2])

    # Save the best_circles to a table.
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.flatfork_asciiname)
    table = Table()
    table['clust_to_try'] = clusters_to_try
    thetapoles, phipoles = np.array(poles_found).T
    table['theta'] = thetapoles
    table['phi'] = phipoles
    writer.write_table("greatcircles_preliminary", "greatcircles", table)

