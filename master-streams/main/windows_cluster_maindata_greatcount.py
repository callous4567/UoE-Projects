import pickle
import numpy as np
from matplotlib import tri
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
from galcentricutils import greatcount
import matplotlib.pyplot as plt

# Load the "mean clustering"
with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(file=f)

# Get the unique clusters in this
clusters_to_try = []
for clust in clustered:
    if clust not in clusters_to_try:
        clusters_to_try.append(clust)



# Set up data
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
table = galcentricutils.angular().get_polar(table)

# Example cluster to try
which_try = 0

# Demarcate which ones to use
included = [True if d == which_try else False for d in clustered]
table_included = table[included]

# Generate the grid of great circles to try for the "first pass"
grid_edge_length = 100
pole_thetas, pole_phis, pole_indices = greatcount().fullgridgen(np.rint(grid_edge_length**2, casting='unsafe'))

# Generate run class
greatfit = galcentricutils.greatfit()
resolution = 100
least_squares = []
for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
    leastsq = greatfit.least_squares(table['theta'], table['phi'], theta_pole, phi_pole, resolution)
    least_squares.append(leastsq)

# Get argmin
best_circle = np.argmin(least_squares)
best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])

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
    plt.show()

# Go forth and plot the greatcircle against the data, in galcentricpolars. Already sorted: don't need "clustering."
graphutils.spec_graph().clust_thetaphi(table=table, clustering=clustered, cluster_id=which_try,
                                       vasiliev=False, savedexdir="testing_thetaphi",
                                       gcc_thetaphis=greatfit.gcc_gen(500, *best_circle))