import math

from astropy.table import Table
from awsimple import sns

import ascii_info
import galcentricutils
import matplotlib as mpl

import graphutils
import hdfutils
import windows_directories
import energistics
import numpy as np
import matplotlib.pyplot as plt
import mpl_scatter_density # adds projection='scatter_density'
from scipy.spatial import cKDTree

"""
panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                             ascii_info.set_raw) """




#energistics.energistics_manual().pot_eval(table)
# https://stackoverflow.com/questions/2369492/generate-a-heatmap-in-matplotlib-using-a-scatter-data-set
def data_coord2view_coord(p, vlen, pmin, pmax):
    dp = pmax - pmin
    dv = (p - pmin) / dp * vlen
    return dv
def nearest_neighbours(xs, ys, reso, n_neighbours):
    im = np.zeros([reso, reso])
    extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]

    xv = data_coord2view_coord(xs, reso, extent[0], extent[1])
    yv = data_coord2view_coord(ys, reso, extent[2], extent[3])
    for x in range(reso):
        for y in range(reso):
            xp = (xv - x)
            yp = (yv - y)

            d = np.sqrt(xp**2 + yp**2)

            im[y][x] = 1 / np.sum(d[np.argpartition(d.ravel(), n_neighbours)[:n_neighbours]])

    return im, extent
def finrun():
    xs = table['circ']
    ys = table['Lz']
    resolution = 400
    fig, ax = plt.subplots(1, 1)
    im, extent = nearest_neighbours(xs, ys, resolution, 16)
    #h = ax.imshow(im, origin='lower', extent=extent, cmap=cm.inferno, aspect='auto')
    h = ax.scatter(xs, ys, s=0.5, cmap=cm.inferno, marker=".")
    fig.colorbar(h, label="Count")
    ax.set_title("Smoothing over %d neighbours" % 4)
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-5e3,5e3)
    plt.show()

# https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
# Alternative method that we will use- more suitable.

# Make the norm object to define the image stretch

from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                             ascii_info.set_raw)

graphutils.spec_graph().sofie(table)

