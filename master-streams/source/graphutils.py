import copy
import os

import astropy
import numpy as np
from astropy.table import Table, unique, vstack

# Our own stuff.
import galcentricutils
import windows_directories

# Matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import matplotlib.colors as mcolors


# Plotly
import plotly.express as px
import plotly.io as plio
plio.renderers.default = 'browser'

# Pandas
import pandas


# Holds all graphing utilities for project. Notes:
"""
Units should be as standard.
xyz in kpc
vxvyvz in kms^-1
lxlylz in kpc kms^-1
angles in degrees
"""

# 2D Graphing Tools
class twod_graph(object):
    def __init__(self):
        self.null = "null"

    # Creates LxLyLxLzLyLz Plots. Axis units of 10^3 kpckms^-1
    def tripL(self, table, maxradius, divisor=1e3):
        table = galcentricutils.angular().get_momentum(table)
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey="row", figsize=(15,7), constrained_layout=True)
        plt.subplots_adjust(wspace=0, hspace=0)

        minor_ticks_x = pltick.MultipleLocator(1)
        minor_ticks_y = pltick.MultipleLocator(1)

        texts = ["Lx - Ly", "Lx - Lz", "Ly - Lz"]
        lims = [-maxradius,maxradius]

        for num,ax in enumerate(axs):
            ax.grid(True, which='major', alpha=0, linestyle='dotted')  # Enable grids on subplot
            ax.grid(True, which='minor', alpha=0, linestyle='dotted')
            ax.set(xlim=lims,
                   ylim=lims)
            ax.set(aspect="equal")
            ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                           top=True)
            ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                           top=True)
            ax.xaxis.set_minor_locator(minor_ticks_x)
            ax.yaxis.set_minor_locator(minor_ticks_y)
            ax.text(s=texts[num],x=lims[0],y=lims[1]+0.1)

        axs[0].scatter(table['Lx']/divisor, table['Ly']/divisor, s=1, marker="x",c="black")
        axs[1].scatter(table['Lx']/divisor,table['Lz']/divisor,s=1,marker="x",c="black")
        axs[2].scatter(table['Ly']/divisor, table['Lz']/divisor, s=1, marker="x",c="black")

        axs[1].set(xlabel="(" + r'$10^3$' + " kpc km/s)")
        axs[0].set(ylabel="(" + r'$10^3$' + " kpc km/s)")

        plt.savefig(windows_directories.imgdir + "\\test.png", dpi=300)
        plt.show()

    # Hammer-Aitoff Projection Plots. Need angles in RADIANS for some god damn ridiculous reason.
    # [-180,180],[-90,90] phi,theta
    def haitoff(self, table):
        # First get latitude and azimuth
        table = galcentricutils.angular().get_latipolar(table)
        phis, thetas = np.deg2rad(table['phi']),np.deg2rad(table['theta'])
        for num,phi in enumerate(phis):
            if phi >= np.pi:
                phis[num] = phi - 2*np.pi
        fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': "aitoff"}, dpi=300)
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
        axs.scatter(phis, thetas, s=0.1, alpha=1, color="black")
        plt.show()

    # 2D Histogram Plot in (l,b)
    def lbtwodhist(self, table, nbins):
        # First get latitude and azimuth
        table = galcentricutils.angular().get_latipolar(table)
        phis, thetas = np.deg2rad(table['phi']),np.deg2rad(table['theta'])
        for num,phi in enumerate(phis):
            if phi >= np.pi:
                phis[num] = phi - 2*np.pi
        fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
        h = axs.hist2d(phis, thetas, bins=[nbins,nbins], cmin=1, cmap = plt.cm.nipy_spectral)
        plt.colorbar(h[3], label="Bin Count")
        axs.grid(False)
        axs.set(xlabel='l',
                ylabel='b')
        plt.show()

    # Regular 2D l,b plot.
    def lbplot(self, table):
        tab_x,tab_y = table['l'],table['b']
        fig, axs = plt.subplots(nrows=1,ncols=1, dpi=300)
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
        axs.scatter(tab_x,tab_y,color="black",s=0.1)
        plt.gca().invert_xaxis()
        axs.set(xlabel="l", ylabel="b")
        plt.show()

    # Regular 2D ra,dec plot.
    def radecplot(self, table):
        tab_x,tab_y = table['ra'],table['dec']
        fig, axs = plt.subplots(nrows=1,ncols=1, dpi=300)
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
        axs.scatter(tab_x,tab_y,color="black",s=0.1)
        plt.gca().invert_xaxis()
        axs.set(xlabel=r'$\alpha$', ylabel=r'$\delta$')
        plt.show()

    # Create a histogram plot for data.
    def hist(self, data, nbins, range):
        if range == False:
            plt.hist(data, bins=nbins, density=True)
            plt.show()
        else:
            plt.hist(data, bins=nbins, density=True, range=range)
            plt.show()

# 3D graphing tools
class threed_graph(object):
    def __init__(self):
        self.null = "null"


    # 3D K-Means Angular Momentum Plot - coloured points. Assumes L exists in table (no resolve.)
    # The "kmeans_L" label is a bit deprecated- we'll roll with it for now, though.
    """
    Same as twod plot for LXYZ. 
    Divisor takes units into /10^3 kpckms^-1 
    """
    def kmeans_L(self, table, savedex, browser, divisor=1e3):
        """
        # Set up figure/etc
        fig = plt.figure(figsize=(15,7), constrained_layout=True)
        ax = fig.add_subplot(projection='3d')
        plt.subplots_adjust(wspace=0, hspace=0)

        # note a list of usable colours.
        colours = ['b', 'g', 'r', 'c', 'm', 'y', 'orange']

        # Set up params for the grid/axis.

        ax.grid(True, which='major', alpha=0, linestyle='dotted')  # Enable grids on subplot
        ax.grid(True, which='minor', alpha=0, linestyle='dotted')

        ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                       top=True)
        ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                       top=True)


        # Get all unique cluster indices for the table
        unique_clusters = unique(table, keys="k_index")['k_index']
        # Get table clips for each unique cluster
        table_by_cluster = table.group_by("k_index")
        tables_for_clust = []
        for cluster in unique_clusters:
            cluster_mask = table_by_cluster.groups.keys["k_index"] == cluster
            masked_table = table_by_cluster.groups[cluster_mask]
            tables_for_clust.append(masked_table)

        # For each cluster, plot in 3D the angular momentum, with colour
        for num, table_for_clust in enumerate(tables_for_clust):
            Lx,Ly,Lz = table_for_clust['Lx'], table_for_clust['Ly'], table_for_clust['Lz']
            ax.scatter(xs=Lx, ys=Ly, zs=Lz, color=colours[num], s=30)

        plt.savefig(windows_directories.imgdir + "\\test3d.png", dpi=300)
        plt.show()

        """

        # Also generate an interactive 3D plot. This is mainly for debug purposes.
        try:
            os.mkdir(windows_directories.imgdir + "\\kmeans_html")
        except:
            pass

        table_pandas = table.to_pandas()
        fig = px.scatter_3d(table_pandas, x='Lx', y='Ly', z='Lz', color='k_index')
        save_loc = windows_directories.imgdir + "\\kmeans_html" + "\\" + savedex
        fig.write_html(save_loc)
        if browser == True:
            fig.show()

    # Basically kmeans_L but spatially instead of in L-space
    def xmeans_L(self, table, savedex, browser, divisor=1e3):
        """
        # Set up figure/etc
        fig = plt.figure(figsize=(15,7), constrained_layout=True)
        ax = fig.add_subplot(projection='3d')
        plt.subplots_adjust(wspace=0, hspace=0)

        # note a list of usable colours.
        colours = ['b', 'g', 'r', 'c', 'm', 'y', 'orange']

        # Set up params for the grid/axis.

        ax.grid(True, which='major', alpha=0, linestyle='dotted')  # Enable grids on subplot
        ax.grid(True, which='minor', alpha=0, linestyle='dotted')

        ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                       top=True)
        ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                       top=True)


        # Get all unique cluster indices for the table
        unique_clusters = unique(table, keys="k_index")['k_index']
        # Get table clips for each unique cluster
        table_by_cluster = table.group_by("k_index")
        tables_for_clust = []
        for cluster in unique_clusters:
            cluster_mask = table_by_cluster.groups.keys["k_index"] == cluster
            masked_table = table_by_cluster.groups[cluster_mask]
            tables_for_clust.append(masked_table)

        # For each cluster, plot in 3D the angular momentum, with colour
        for num, table_for_clust in enumerate(tables_for_clust):
            Lx,Ly,Lz = table_for_clust['Lx'], table_for_clust['Ly'], table_for_clust['Lz']
            ax.scatter(xs=Lx, ys=Ly, zs=Lz, color=colours[num], s=30)

        plt.savefig(windows_directories.imgdir + "\\test3d.png", dpi=300)
        plt.show()

        """

        # Also generate an interactive 3D plot. This is mainly for debug purposes.
        try:
            os.mkdir(windows_directories.imgdir + "\\xmeans_html")
        except:
            pass

        table_pandas = table.to_pandas()
        fig = px.scatter_3d(table_pandas, x='x', y='y', z='z', color='k_index', title=savedex)
        save_loc = windows_directories.imgdir + "\\xmeans_html" + "\\" + savedex
        fig.write_html(save_loc)
        if browser == True:
            fig.show()

    # Given a set of theta,phi plot on the unit sphere in 3D for browser. "Colours" are optional.
    def unitsphere(self, thetas, phis, colours):
        pos = lambda theta, phi: np.array([np.sin(theta)*np.cos(phi),
                                           np.sin(theta)*np.sin(phi),
                                           np.cos(theta)])

        thetaphis = list(zip(thetas,phis))
        coords = []
        for thetaphi in thetaphis:
            coords.append(pos(thetaphi[0],thetaphi[1]))
        posses = np.array(coords)
        posses = posses.T

        postable = Table()
        postable['x'],postable['y'],postable['z'] = posses

        if type(colours) != bool:
            postable['colour'] = colours
        else:
            postable['colour'] = ["400" for d in postable['z']]
        table_pandas = postable.to_pandas()
        fig = px.scatter_3d(table_pandas, x='x', y='y', z='z', color="colour")
        fig.show()