import pickle
import numpy as np
from astropy.table import Table
from matplotlib import tri
import ascii_info_new
import galcentricutils_new
import graphutils_new
import hdfutils
import windows_directories_new
from energistics_new import energistics_new, energistics_new_manual, orbigistics
from galcentricutils_new import greatcount
import matplotlib.pyplot as plt

# Load the post-initial-montecarlo-clustering (via memberpercent.py.)
"""
This has the benefit of having weeded-out most of the "non-clusterings" via monte-carlo-
most of them will have vanished by this point and been absorbed into other clusters, if they're not the "significant one."
"""
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
percent_table = writer.read_table(ascii_info_new.fullgroup, "percent_table")
clustered = percent_table['probable_clust']

#with open(windows_directories_new.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
#    clustered = pickle.load(file=f)

# Get the unique clusters in this
clusters_to_try = []
for clust in clustered:
    if clust not in clusters_to_try:
        clusters_to_try.append(clust)

# Set up data
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.fullset)
table = galcentricutils_new.angular().get_polar(table)

# Generate the grid of great circles to try for the "first pass"
grid_edge_length = 400 # number of GCC poles to generate
resolution = 800 # number of datapoints to generate in GCC
pole_thetas, pole_phis, pole_indices = greatcount().fullgridgen(np.rint(grid_edge_length**2, casting='unsafe'))
pole_thetas, pole_phis = np.degrees(pole_thetas), np.degrees(pole_phis)

# Generate run class
greatfit = galcentricutils_new.greatfit()

# Do this method if you're dealing with "real_dist=False" and just least-squaring.
"""
real_dist will specify the distance for which to maximize enclosed points to the great circle: see the implementation.
"""
real_dist = 5

# Poles Found
poles_found = []


# For each individual clustering we have, try to fit greatcircles and produce plots
for clust in clusters_to_try:
    # Little catching specification for specific clusters (specifically, the orphan.)
    #if clust == 1: # Orphan
    #    real_dist_local = 2
    #else:
    real_dist_local = real_dist

    # Example cluster to try
    which_try = clust

    # Demarcate which ones to use
    included = [True if d == which_try else False for d in clustered]
    table_included = table[included]

    if real_dist_local == False:
        least_squares = []
        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole, resolution, real_dist_local)
            least_squares.append(leastsq)
        best_circle = np.argmin(least_squares)
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])
        poles_found.append(best_circle)
    # Else use the radius you set.
    else:
        least_squares = []
        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,resolution, real_dist_local)
            least_squares.append(leastsq)
        best_circle = np.argmax(least_squares) # argmax, to mamimize the membership of these cluster members (using a real_dist- see greatfit.least_squares())
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])
        poles_found.append(best_circle)
    """
    # Create contour plot
    contour = False
    if contour == True:
        ngrid = 100
        xi = np.linspace(np.min(pole_thetas), np.max(pole_thetas), ngrid)
        yi = np.linspace(np.min(pole_phis), np.max(pole_phis), ngrid)
        triang = tri.Triangulation(pole_thetas, pole_phis)
        interpolator = tri.LinearTriInterpolator(triang, least_squares)
        Xi, Yi = np.meshgrid(xi, yi)
        zi = interpolator(Xi, Yi)
    
        # Get fig
        fig, ax1 = plt.subplots(nrows=1, ncols=1)
        ax1.contour(xi, yi, zi, levels=14, linewidths=0.5, colors='k')
        cntr1 = ax1.contourf(xi, yi, zi, levels=14, cmap="RdBu_r")
        fig.colorbar(cntr1, ax=ax1)
        plt.show() """

    # Go forth and plot the greatcircle against the data, in ra/dec.
    graphutils_new.spec_graph().clust_thetaphi(table=table_included, clustering=clustered[included], cluster_id=which_try,
                                           vasiliev=False, savedexdir="greatcount_" + str(clust) + "_thetaphi_greatcircle",
                                           gcc_thetaphis=greatfit.gcc_gen(1000, *best_circle))

# Save the best_circles to a table.
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
table = Table()
table['clust_to_try'] = clusters_to_try
thetapoles, phipoles = np.array(poles_found).T
table['theta'] = thetapoles
table['phi'] = phipoles
writer.write_table("greatcircles_preliminary", "greatcircles", table)

