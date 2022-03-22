import copy
import itertools
import os
import time

import astropy
import numpy as np
from scipy.stats import gaussian_kde
from astropy.table import Table, unique, vstack
from astropy import units as u

# Our own stuff.
from galpy.util import bovy_coords
from matplotlib.patches import Patch

import galcentricutils
import hdfutils
import windows_directories

# Matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import matplotlib.colors as mcolors
from matplotlib import cm, rc, colors

# Plotly
import plotly.express as px
import plotly.io as plio



plio.renderers.default = 'browser'
import plotly.graph_objects as go

# Pandas
import pandas

# Pillow
import PIL
from PIL import Image

# Misc Astropy
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import units as u

# Holds all graphing utilities for project. Notes:
"""
Units should be as standard.
xyz in kpc
vxvyvz in kms^-1
lxlylz in kpc kms^-1
angles in degrees
"""

# Enable TeX
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# 2D Graphing Tools prototyping
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

    # Latipolar plot without Hammer Projection, with a great circle [theta,phi] plotted.
    # Theta,phi are in regular galactocentric polar
    def latipolar(self, table, theta_polar, phi):
        # First get latitude and azimuth
        table = galcentricutils.angular().get_latipolar(table)
        phis, thetas = np.deg2rad(table['phi']),np.deg2rad(table['theta'])
        fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
        phis, thetas = np.rad2deg(phis), np.rad2deg(thetas)
        axs.scatter(phis, thetas, s=0.1, alpha=1, color="black")

        # Add the greatcircle you've opted to display
        gcc = galcentricutils.greatcount().gcc_gen(10000, theta_polar, phi) # thetas, phis
        plt.scatter(gcc[1],90 - gcc[0], color="blue", s=1)

        plt.show()

    # 2D Histogram Plot in (theta,phi): calculates polar for the table.
    def thetaphi_twodhist(self, table, nbins, savedex, latipolar):
        # First get latitude and azimuth
        table = galcentricutils.angular().get_polar(table)
        phis, thetas = table['phi'],table['theta']
        #for num,phi in enumerate(phis):
        #    if phi >= np.pi:
        #        phis[num] = phi - 2*np.pi
        if latipolar == True:
            thetas = 90 - thetas
        fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
        h = axs.hist2d(phis, thetas, bins=[nbins,nbins], cmin=1, cmap = plt.cm.nipy_spectral)
        plt.colorbar(h[3], label="Bin Count")
        axs.grid(False)
        axs.set(xlabel=r'$\phi$',
                ylabel=r'$\theta$')
        if latipolar == True:
            axs.set(ylabel=r'$\theta$' + " in Latipolar")
        try:
            os.mkdir(windows_directories.imgdir + "\\thetaphi_twodhist")
        except:
            pass

        plt.savefig(windows_directories.imgdir + "\\thetaphi_twodhist" + "\\" + savedex + ".png")
        plt.show()

    # Fancy l,b plot
    def lbplot(self, table, savepath, negpi=False):

        # Get axes and data
        plt.close()
        tab_x,tab_y = table['l'],table['b']
        if negpi == True:
            tab_x = [d - 360 if d > 180 else d for d in tab_x]
        fig, axs = plt.subplots(nrows=1,ncols=1, dpi=300)
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="pink")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="pink")
        axs.scatter(tab_x,tab_y,color="red",s=0.5)
        plt.gca().invert_xaxis()
        axs.set(xlabel="l / deg", ylabel="b / deg")
        plt.savefig(savepath, dpi=300)
        plt.close()

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

    # Plot nclust vs. n for a monte-carloed clustering set.
    # Also subplot the PDF for the clustering.
    def nclust_n(self, data, group):
        # The n_clust/clustering plot.
        x = np.arange(0, len(data), 1)
        y = data
        fig, axs = plt.subplots(nrows=1,ncols=2,sharey='row', gridspec_kw={'width_ratios': [3,1]})
        plt.subplots_adjust(wspace=0.1, hspace=0)
        axs[0].scatter(x, y, marker='x',color='black',s=10)
        axs[0].set(xlabel="monte-carlo clustering/index",
               ylabel="n / clusters")

        # Generate a PDF, too.
        axs[1].hist(y, bins=60, density=True, orientation='horizontal', color='black')
        axs[1].set(xlabel=r'$\rho$')
        #plt.suptitle("For " + group.replace("_","\_"))
        plt.savefig(windows_directories.imgdir + "\\graph_nclust.png", dpi=300)
        plt.show()

    # Create a very fancy histogram plot for the data! (PDF.)
    def hist_fancy(self, data, nbins, xlabel, ylabel, savepath=None, dencurve=True):

        # Figure
        fig = plt.figure()
        # Set up the grid
        plt.grid(which='major', color='pink')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # Without smooth curve to plot
        if dencurve == False:

            # With/without range
            plt.hist(data, bins=nbins, density=True, color='lightblue')

        # With smooth curve to plot over the histogram
        else:
            # Do a hist
            plt.hist(data, bins=nbins, density=True, color='lightblue')

            # Generate smooth plot to pdf
            kde = gaussian_kde(data)
            xvalues = np.linspace(np.min(data), np.max(data), 1000)
            yvalues = kde.evaluate(xvalues)
            plt.plot(xvalues, yvalues, lw=0.5, color='red')


        # Save the figure or plot it
        if savepath == None:
            plt.show()
        else:
            plt.savefig(savepath + ".png", dpi=300)

# 3D graphing tools prototyping
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
        if savedex != False:
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

    # Except it takes arrays. Savedexdir here should be the directory\\savename
    def xmeans_L_array(self, array, clusterdata, savedexdir, browser, outliers=False, omit=None):
        if outliers == False:
            outlier_index = np.where(clusterdata == -1)[0]
            array = np.delete(array, outlier_index, axis=0)
            clusterdata = np.delete(clusterdata, outlier_index, axis=0)

        # Need arrays. Make sure.
        if type(array) != "numpy.ndarray":
            array = np.array(array)
        if type(clusterdata) != "numpy.ndarray":
            clusterdata = np.array(clusterdata)

        # If omit isn't none, then remove omitted clusters
        if omit != None:
            for i in omit:
                truefalse = [False if d == i else True for d in clusterdata]
                clusterdata = clusterdata[truefalse]
                array = array[truefalse]

        x,y,z = array.T

        fig = px.scatter_3d(x=x, y=y, z=z, color=clusterdata)
        if savedexdir != False :
            save_loc = windows_directories.imgdir + "\\" + savedexdir + ".html"
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

    # Except it takes arrays. Savedexdir here should be the directory\\savename
    def kmeans_L_array(self, array, clusterdata, savedexdir, browser, outliers=False):
        # Need arrays. Make sure.
        if type(array) != "numpy.ndarray":
            array = np.array(array)
        if type(clusterdata) != "numpy.ndarray":
            clusterdata = np.array(clusterdata)

        # Verify if array is only Lx Ly Lz or if it also energy, and if so, then clip to just angular momentum.
        if len(array[0]) > 3:
            array = array.T
            array = array[0:3]
            array = array.T
        else:
            pass

        # Decide if you want to keep or remove the noise.
        if outliers == False:
            outlier_index = np.where(clusterdata == -1)[0]
            array = np.delete(array, outlier_index, axis=0)
            clusterdata = np.delete(clusterdata, outlier_index, axis=0)
        else:
            # Limit to the viewing region
            array_clip_indices = galcentricutils.cluster3d().getCuboid(array)
            array, clusterdata = array[array_clip_indices], clusterdata[array_clip_indices]

        x,y,z = array.T

        fig = px.scatter_3d(x=x, y=y, z=z, color=clusterdata)
        if savedexdir != False :
            save_loc = windows_directories.imgdir + "\\" + savedexdir + ".html"
            fig.write_html(save_loc)
        if browser == True:
            fig.show()

    # Given two clusterings [data1,data2] and label sets [labels1,labels2] plot both
    # Shifts the 2nd labelling by 16 for colours.
    def kmeans_L_array_pair(self, arrays, clusterdatas, savedexdir=False, browser=False, outliers=False):
        # If outliers=False, then go and remove all (-1) elements from both arrays.
        if outliers == False:
            for num, clusterdata in enumerate(clusterdatas):
                outlier_index = np.where(clusterdata == -1)[0]
                array = arrays[num]
                arrays[num] = np.delete(array, outlier_index, axis=0)
                clusterdatas[num] = np.delete(clusterdata, outlier_index, axis=0)

        # Need arrays. Make sure.
        for num, array in enumerate(arrays):
            if type(array) != "numpy.ndarray":
                arrays[num] = np.array(array)
            if type(clusterdatas[num]) != "numpy.ndarray":
                clusterdatas[num] = np.array(clusterdatas[num])

        # Limit to the viewing region
        for num, array in enumerate(arrays):
            array_clip_indices = galcentricutils.cluster3d().getCuboid(array)
            arrays[num], clusterdatas[num] = array[array_clip_indices], clusterdatas[num][array_clip_indices]

        # Apply colour shift to visualize cluster pairing
        clusterdatas[1] = [d + 16 for d in clusterdatas[1]]
        # Do the etc etc graphing
        array = np.concatenate(arrays)
        clusterdata = np.concatenate(clusterdatas)
        x,y,z = array.T
        fig = px.scatter_3d(x=x, y=y, z=z, color=clusterdata)
        if savedexdir != False :
            save_loc = windows_directories.imgdir + "\\" + savedexdir + ".html"
            fig.write_html(save_loc)
        if browser == True:
            fig.show()

    # Plots data alongside 3D gaussian convex hulls for provided gaussians. Data as vertical np array
    # [[v1],
    #  [v2],
    #  ....]
    # npoints for convex hull of the multinormals
    # Does not support colouring datasets individually.
    def kmeans_L_multinormal(self, vec_L, mus, covtrices, npoints):
        # Define rng
        rng = np.random.default_rng()
        # Tracelist
        gotracelist = []
        # Set up trace for the actual data
        vecx,vecy,vecz = vec_L.T
        vectrace = go.Scatter3d(x=vecx,
                                y=vecy,
                                z=vecz,
                                name="Data",
                                marker=dict(color='rgb(34,163,192)'),
                                mode='markers')
        gotracelist.append(vectrace)

        # Set up traces for each mu/covtrix (we're not trying to be fancy here- just somewhat dirty.
        for num, mu, cov in zip(range(len(mus)),mus,covtrices):
            points = rng.multivariate_normal(mean=mu,cov=cov,size=npoints)
            points = points.T
            mvtrace = go.Mesh3d(x=points[0],y=points[1],z=points[2],
                                alphahull=0,opacity=0.25,
                                name=("MV {}").format(str(num)))
            gotracelist.append(mvtrace)

        fig = go.Figure(data=gotracelist)
        fig.show()

    # Visualize generated multivariate normals from noisy_data.
    # n_mesh = number of points per mesh (for the mucov multivariate normal.)
    # n_points = a list or tuple with the number of RNG per normal
    def kmeans_L_multinormal_generated(self,mucovs,n_points,n_mesh):
        # Mus and Covs as in generated multinormals
        mus, covtrices = mucovs

        # Define rng
        rng = np.random.default_rng()

        # Tracelist
        gotracelist = []

        # Set up traces for each mu/covtrix (we're not trying to be fancy here- just somewhat dirty.
        for num, mu, cov, npoint in zip(range(len(mus)), mus, covtrices, n_points):
            # Generate convex hull
            points = rng.multivariate_normal(mean=mu,cov=cov,size=n_mesh)
            points = points.T
            meshtrace = go.Mesh3d(x=points[0],y=points[1],z=points[2],
                                  alphahull=0,opacity=0.25,
                                  name=("Mesh {}").format(str(num)),
                                  color=int(num))
            points = points[0:npoint]
            scattrace = go.Scatter3d(x=points[0],y=points[1],z=points[2],
                                     name=("Scatter {}").format(str(num)),
                                     mode='markers')
            gotracelist.append(meshtrace)
            gotracelist.append(scattrace)

        fig = go.Figure(data=gotracelist)
        fig.show()

# More specific graphing functions that may make use of/are based on the above, or be novel.
class spec_graph(object):
    def __init__(self):
        self.null = "null"

    """
    radec plot. 
    Select a cluster_id, otherwise all are plotted. 
    Set vasiliev to True to display his Sgr_snapshot (with LMC.) 
    optionally supply gcc_radecs = [[ras],[decs]] with greatcircle radecs
    """
    def clust_radec(self, table, clustering, cluster_id=False, vasiliev=False, savedexdir=None, gcc_radecs=None, lb=False):
        if cluster_id is False:
            # Get the n0 of clusters
            nclust = np.max(clustering) + 1

            # Also grab the order of the clusters, too, based on which is largest/smallest.

            # ICRS-up the thing and do the plot...
            fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)

            # first do Vasilievs, though.
            if vasiliev==True:
                simdir = windows_directories.datadir + "\\vasiliev\\Sgr_snapshot"
                simfile = "simdata.hdf5"
                writer = hdfutils.hdf5_writer(simdir, simfile)
                sim_table = writer.read_table("Sgr_snapshot", "astrotable")
                if lb == False:
                    tab_x, tab_y = sim_table['ra'], sim_table['dec']
                    axs.scatter(tab_x, tab_y, color="gray", s=0.01, marker=".")
                elif lb == True:
                    l, b = bovy_coords.radec_to_lb(sim_table['ra'], sim_table['dec'], degree=True).T
                    l = [d - 360 if d > 180 else d for d in l]
                    axs.scatter(l, b, color="gray", s=0.01, marker=".")

            # Sort table by clustering
            table = galcentricutils.galconversion().nowrite_GAL_to_ICRS(table)
            table['cluster'] = clustering
            table_by_clust = table.group_by("cluster")

            # Miscellaneous thing to ensure non-cyclical rainbow colourmap.
            cm = plt.get_cmap('gist_rainbow')
            axs.set_prop_cycle('color', [cm(1. * i / nclust) for i in range(nclust)])

            # Define an arbitrary function to set the scale based on the number of stars in the clustering.
            # This is just to ensure that the larger clusters don't just block out the smaller ones.
            scale = lambda n: 0.1 + 5*((30/n)**0.8)

            # Save the masked tables
            masked_tables = []
            masked_tables_lengths = []

            # Get masked tables for each cluster, alongside the number of stars in each
            for cluster_id in range(nclust):
                clustmask = table_by_clust.groups.keys["cluster"] == cluster_id
                masked_table = table_by_clust.groups[clustmask]
                masked_tables.append(masked_table)
                masked_tables_lengths.append(len(masked_table))

            # Sort the tables and table scales (just to ensure we don't hide away small ones under big ones)
            sorted_tables = [x for (y, x) in sorted(zip(masked_tables_lengths, masked_tables), key=lambda pair: pair[0])]
            sorted_tables.reverse()

            # Now, plot for the sorted tables
            for sorted_table in sorted_tables:
                if lb == False:
                    plt.gca().invert_xaxis()
                    tab_x, tab_y = sorted_table['ra'], sorted_table['dec']
                    axs.scatter(tab_x, tab_y, s=scale(len(sorted_table)), marker="s")
                elif lb == True:
                    l, b = bovy_coords.radec_to_lb(sorted_table['ra'], sorted_table['dec'], degree=True).T
                    l = [d - 360 if d > 180 else d for d in l]
                    axs.scatter(l, b, color="gray", s=0.01, marker=".")
            plt.gca().invert_xaxis()

            # Misc axis things
            axs.set_facecolor("k")
            axs.grid(True, which='major', alpha=1, linewidth=0.25)
            axs.grid(True, which='minor', alpha=1, linewidth=0.25)
            axs.grid(color="white")

            if lb == False:
                axs.set(xlabel=r'$\alpha$' + " / deg", ylabel=r'$\delta$' + " / deg")
            elif lb == True:
                axs.set(xlabel=" l / deg", ylabel="b / deg")

            # Create save directory/save
            if savedexdir != None:
                try:
                    os.mkdir(windows_directories.imgdir + "\\" + "vasiliev")
                except:
                    pass
                plt.savefig(windows_directories.imgdir + "\\" + "vasiliev"+ "\\" + savedexdir + ".png", dpi=300)
            # Show
            #plt.show()

        else:
            # Grab the ra/dec list for the cluster_id
            table['cluster'] = clustering
            table_by_clust = table.group_by("cluster")
            clustmask = table_by_clust.groups.keys["cluster"] == cluster_id
            masked_table = table_by_clust.groups[clustmask]

            # ICRS-up the thing and do the plot...
            fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
            # first do Vasilievs, though.
            legend_elements = []
            if vasiliev==True:
                legend_elements.append(Patch(edgecolor='white',facecolor='gray',label="Vasiliev Model"))
                simdir = windows_directories.datadir + "\\vasiliev\\Sgr_snapshot"
                simfile = "simdata.hdf5"
                writer = hdfutils.hdf5_writer(simdir, simfile)
                sim_table = writer.read_table("Sgr_snapshot", "astrotable")
                if lb == False:
                    tab_x, tab_y = sim_table['ra'], sim_table['dec']
                    axs.scatter(tab_x, tab_y, color="gray", s=0.01, marker=".")
                elif lb == True:
                    l, b = bovy_coords.lb_to_radec(sim_table['ra'], sim_table['dec'], degree=True).T
                    l = [d - 360 if d > 180 else d for d in l]
                    axs.scatter(l, b, color="gray", s=0.01, marker=".")
            # Then our clustering...
            masked_table = galcentricutils.galconversion().nowrite_GAL_to_ICRS(masked_table)

            if lb == False:
                tab_x, tab_y = masked_table['ra'], masked_table['dec']
                axs.scatter(tab_x, tab_y, s=5, marker="s", color='red')
            elif lb == True:
                l, b = bovy_coords.radec_to_lb(masked_table['ra'], masked_table['dec'], degree=True).T
                l = [d - 360 if d > 180 else d for d in l]
                axs.scatter(l, b, color="red", s=1, marker="o")

            axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
            axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
            legend_elements.append(Patch(edgecolor='white',
                                         facecolor='red',
                                         label="Stars"))
            # We can also add a great circle (though this necessitates a legend.)
            if gcc_radecs != None:
                if lb == False:
                    gccras, gccdecs = gcc_radecs
                    axs.scatter(gccras, gccdecs, color="blue", s=0.1)
                elif lb == True:
                    gccras, gccdecs = gcc_radecs
                    l, b = bovy_coords.radec_to_lb(gccras, gccdecs, degree=True).T
                    axs.scatter(l, b, color="blue", s=0.1)
                legend_elements.append(Patch(edgecolor='white',facecolor='blue',label='GCC Fit'))
            plt.legend(loc='upper right', handles=legend_elements)
            axs.set_facecolor("k")
            if lb == False:
                axs.set(xlabel=r'$\alpha$' + " / deg", ylabel=r'$\delta$' + " / deg")
            elif lb == True:
                axs.set(xlabel=" l / deg", ylabel="b / deg")
            axs.grid(color="white")
            plt.gca().invert_xaxis()

            # Create save directory/save
            if savedexdir != None:
                try:
                    os.mkdir(windows_directories.imgdir + "\\" + "vasiliev")
                except:
                    pass
                plt.savefig(windows_directories.imgdir + "\\" + "vasiliev"+ "\\" + savedexdir + ".png", dpi=300)
            plt.close()
            #plt.show()

    # Colour-coded map like in "The Global Dynamical Atlas of the Milky Way Mergers:
    # Constraints from Gaia EDR3–based Orbits of Globular Clusters, Stellar Streams, and Satellite Galaxies"
    # Coded by proper motion
    def clust_lb_pms(self, table, clustering, clusts_to_plot, savepath):

        # Create axis
        fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(10,4))

        def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=1000):
            new_cmap = colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)))
            return new_cmap

        # Set up colormap
        cmap = cm.get_cmap('hsv_r')
        cmap = truncate_colormap(cmap, 0.3, 1)
        norm = colors.Normalize(vmin=-15, vmax=15)

        # Plot
        for clust in clusts_to_plot:

            retain = [True if d == clust else False for d in clustering]
            data_to_fit = table[retain]
            l, b = data_to_fit['l'], data_to_fit['b']
            l = [d - 360 if d > 180 else d for d in l]
            dmu_b = data_to_fit['dmu_b']
            plt.scatter(x=l, y=b, s=5, c=dmu_b, cmap=cmap, norm=norm)

        # Grid-up
        # plt.legend(handles=legend_elements, loc='upper right')
        axs.grid(which='major', color='pink')
        axs.set_facecolor("k")
        plt.gca().invert_xaxis()
        plt.colorbar(label=r'$\mu_\textrm{b}$' + " / mas" + r'$\cdot\textrm{yr}^{-1}$')
        if savepath != None:
            try:
                plt.savefig(savepath, dpi=300)
            except:
                plt.savefig(savepath)

        plt.show(dpi=200)

    # The above but with colourmaps
    def clust_lb_colors(self, table, clustering, clusts_to_plot, savepath):

        # Create axis
        fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(7,4))
        axs.set(xlim=[-180,180])
        axs.set(ylim=[-90,90])

        c = ["red",
             "darkorange",
             "lime",
             "cyan",
             "blue",
             "fuchsia",
             "darkviolet",
             "dodgerblue",
             "chocolate"]
        if len(c) < len(clusts_to_plot):
            randoms = np.random.random((len(clusts_to_plot) - len(c), 3))
            for i in randoms:
                c.append(i)

        # Legend elements
        legend_elements = []

        # Plot
        for num,clust in enumerate(clusts_to_plot):

            retain = [True if d == clust else False for d in clustering]
            data_to_fit = table[retain]
            l, b = data_to_fit['l'], data_to_fit['b']
            l = [d - 360 if d > 180 else d for d in l]
            plt.scatter(x=l, y=b, s=5, c=c[num])
            legend_elements.append(Patch(edgecolor='black',
                                         facecolor=c[num],
                                         label=str(clust)))


        # Grid-up
        # plt.legend(handles=legend_elements, loc='upper right')
        axs.grid(which='major', color='pink')
        axs.set_facecolor("k")
        plt.gca().invert_xaxis()
        plt.legend(handles=legend_elements, loc='lower left')
        if savepath != None:
            try:
                plt.savefig(savepath, dpi=300)
            except:
                plt.savefig(savepath)

        plt.show(dpi=200)

    # Create a plot with orbit integrals for a table and display in l/b space
    def lb_orbits(self, table, time_in_years, xlim, ylim, savepath=None, line=False, points=1000):

        # Set axes
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))

        # Make it pretty
        ax.set(xlabel="l / deg",
               ylabel="b / deg",
               xlim=xlim,
               ylim=ylim)

        # Generate all the plots for each orbit
        import energistics
        orbigist = energistics.orbigistics()
        orbits = orbigist.orbits(table)
        inttime = np.linspace(0, time_in_years, points)*u.yr
        for orbit in orbits:
            forward = copy.copy(orbit)
            backward = copy.copy(orbit)
            forward.integrate(inttime, orbigist.pot)
            forward = forward.getOrbit()
            R, vR, vT, z, vz, phi = forward.T
            R *= orbigist.rovo[0]
            vR *= orbigist.rovo[1]
            vT *= orbigist.rovo[1]
            z *= orbigist.rovo[0]
            vz *= orbigist.rovo[1]
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
            l = [d - 360 if d > 180 else d for d in l]
            if line == True:
                ax.plot(l,b,lw=1,color='red', zorder=1)
            else:
                ax.scatter(l, b, s=0.01, color='red', zorder=1)

            backward.integrate(-inttime, orbigist.pot)
            backward = backward.getOrbit()
            R, vR, vT, z, vz, phi = backward.T
            R *= orbigist.rovo[0]
            vR *= orbigist.rovo[1]
            vT *= orbigist.rovo[1]
            z *= orbigist.rovo[0]
            vz *= orbigist.rovo[1]
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
            l = [d - 360 if d > 180 else d for d in l]
            if line == True:
                ax.plot(l,b,lw=1,color='blue', zorder=1)
            else:
                ax.scatter(l, b, s=0.01, color='blue', zorder=1)

        # Plot the base
        ax.scatter([d - 360 if d > 180 else d for d in table['l']], table['b'], s=15, color='white', marker='x',
                   zorder=2)

        ax.set_facecolor("k")
        plt.gca().invert_xaxis()

        if savepath != None:
            try:
                plt.savefig(savepath, dpi=300)
            except:
                plt.savefig(savepath)

        plt.show(dpi=200)

    # above but for theta/phi
    def clust_thetaphi(self, table, clustering, cluster_id=False, vasiliev=False, savedexdir=None, gcc_thetaphis=None):
        if cluster_id is False:
            # Get the n0 of clusters
            nclust = np.max(clustering) + 1
            fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)

            # first do Vasilievs, though.
            if vasiliev==True:
                simdir = windows_directories.datadir + "\\vasiliev\\Sgr_snapshot"
                simfile = "simdata.hdf5"
                writer = hdfutils.hdf5_writer(simdir, simfile)
                sim_table = writer.read_table("Sgr_snapshot", "astrotable")
                # Convert sim_table to galcentric
                converter = galcentricutils.galconversion()
                converter.solinfo_grab(windows_directories.sourcedir, "solar_info.dat")
                converter.solgal_set()
                sim_table = converter.nowrite_GAL_to_GALCENT(sim_table)
                # Get polar coordinates for vasiliev sim
                sim_table = galcentricutils.angular().get_polar(sim_table)
                tab_x, tab_y = sim_table['theta'], sim_table['phi']
                axs.scatter(tab_x, tab_y, color="gray", s=0.01, marker=".")

            # Sort table by clustering (not needed.)
            table['cluster'] = clustering
            table_by_clust = table.group_by("cluster")

            # Miscellaneous thing to ensure non-cyclical rainbow colourmap.
            cm = plt.get_cmap('gist_rainbow')
            axs.set_prop_cycle('color', [cm(1. * i / nclust) for i in range(nclust)])

            # Define an arbitrary function to set the scale based on the number of stars in the clustering.
            # This is just to ensure that the larger clusters don't just block out the smaller ones.
            scale = lambda n: 0.1 + 5*((30/n)**0.8)

            # Save the masked tables
            masked_tables = []
            masked_tables_lengths = []

            # Get masked tables for each cluster, alongside the number of stars in each
            for cluster_id in range(nclust):
                clustmask = table_by_clust.groups.keys["cluster"] == cluster_id
                masked_table = table_by_clust.groups[clustmask]
                masked_tables.append(masked_table)
                masked_tables_lengths.append(len(masked_table))

            # Sort the tables and table scales (just to ensure we don't hide away small ones under big ones)
            sorted_tables = [x for (y, x) in sorted(zip(masked_tables_lengths, masked_tables), key=lambda pair: pair[0])]
            sorted_tables.reverse()

            # Now, plot for the sorted tables
            for sorted_table in sorted_tables:
                tab_x, tab_y = sorted_table['theta'], sorted_table['phi']
                axs.scatter(tab_x, tab_y, s=scale(len(sorted_table)), marker="s")

            # Misc axis things
            axs.set_facecolor("k")
            plt.gca().invert_xaxis()
            axs.grid(True, which='major', alpha=1, linewidth=0.25)
            axs.grid(True, which='minor', alpha=1, linewidth=0.25)
            axs.set(xlabel=r'$\theta$',
                    ylabel=r'$\phi$')
            axs.grid(color="white")

            # Create save directory/save
            if savedexdir != None:
                try:
                    os.mkdir(windows_directories.imgdir + "\\" + "vasiliev")
                except:
                    pass
                plt.savefig(windows_directories.imgdir + "\\" + "vasiliev"+ "\\" + savedexdir + "_thetaphi.png", dpi=300)
            # Show
            #plt.show()

        else:
            # Grab the ra/dec list for the cluster_id
            table['cluster'] = clustering
            table_by_clust = table.group_by("cluster")
            clustmask = table_by_clust.groups.keys["cluster"] == cluster_id
            masked_table = table_by_clust.groups[clustmask]

            # ICRS-up the thing and do the plot...
            fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
            # first do Vasilievs, though.
            legend_elements = []
            if vasiliev==True:
                simdir = windows_directories.datadir + "\\vasiliev\\Sgr_snapshot"
                simfile = "simdata.hdf5"
                writer = hdfutils.hdf5_writer(simdir, simfile)
                sim_table = writer.read_table("Sgr_snapshot", "astrotable")
                # Convert sim_table to galcentric
                converter = galcentricutils.galconversion()
                converter.solinfo_grab(windows_directories.sourcedir, "solar_info.dat")
                converter.solgal_set()
                sim_table = converter.nowrite_GAL_to_GALCENT(sim_table)
                # Get polar coordinates for vasiliev sim
                sim_table = galcentricutils.angular().get_polar(sim_table)
                tab_x, tab_y = sim_table['theta'], sim_table['phi']
                axs.scatter(tab_x, tab_y, color="gray", s=0.01, marker=".")
                legend_elements.append(Patch(edgecolor='white',facecolor='gray',label="Vasiliev Model"))
            # Then our clustering...
            tab_x, tab_y = masked_table['theta'], masked_table['phi']
            axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
            axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
            axs.scatter(tab_x, tab_y, color="red", s=5, marker="s")
            legend_elements.append(Patch(edgecolor='white',
                                         facecolor='red',
                                         label="Stars"))
            # We can also add a great circle (though this necessitates a legend.)
            if gcc_thetaphis != None:
                gccthetas, gccphis = gcc_thetaphis
                axs.scatter(gccthetas, gccphis, color="blue", s=0.1)
                legend_elements.append(Patch(edgecolor='white',facecolor='blue',label='GCC Fit'))

            plt.legend(loc='upper right', handles=legend_elements)
            axs.set_facecolor("k")
            plt.gca().invert_xaxis()
            axs.grid(color="white")
            axs.set(xlabel=r'$\theta$',
                    ylabel=r'$\phi$')
            # Create save directory/save
            if savedexdir != None:
                try:
                    os.mkdir(windows_directories.imgdir + "\\" + "vasiliev")
                except:
                    pass
                plt.savefig(windows_directories.imgdir + "\\" + "vasiliev"+ "\\" + savedexdir + "_thetaphi.png", dpi=300)
            #plt.show()

    # above but for arrays [v1, v2, v3] and will plot ALL points alongside the greatcircle theta/phi. Debug purposes.
    def array_thetaphi_debug(self, array, theta, phi, clusters, whichclust):
        # Get angular
        ang = galcentricutils.angular()
        # Get polars
        polars = [ang.vec_polar(d) for d in array]
        rads,thets,phis = np.array(polars).T
        # Sort plots
        fig, axs = plt.subplots(nrows=1,ncols=1)
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
        # Sort colours/points for the cluster of interest.
        is_cluster = [True if d == whichclust else False for d in clusters]
        is_not = [True if d != whichclust else False for d in clusters]
        # The non members
        axs.scatter(thets[is_not], phis[is_not], color="red", s=5, marker="s")
        # Do the plots for the cluster members
        axs.scatter(thets[is_cluster], phis[is_cluster], color="green", s=3, marker="s")
        legend_elements = []
        legend_elements.append(Patch(edgecolor='white',
                                     facecolor='red',
                                     label="Not Cluster"))
        legend_elements.append(Patch(edgecolor='white',
                                     facecolor='green',
                                     label="In Cluster"))
        # We can also add a great circle (though this necessitates a legend.)
        gccthetas, gccphis = galcentricutils.greatfit().gcc_gen(1000, theta, phi)
        axs.scatter(gccthetas, gccphis, color="blue", s=0.1)
        legend_elements.append(Patch(edgecolor='white',facecolor='blue',label='GCC Fit'))
        plt.legend(loc='upper right', handles=legend_elements)
        axs.set_facecolor("k")
        plt.gca().invert_xaxis()
        axs.grid(color="white")
        axs.set(xlabel=r'$\theta$',
                ylabel=r'$\phi$')
        plt.title(str(whichclust))
        plt.show()

    # Misc func for handling image tiling.
    # https://stackoverflow.com/questions/37921295/python-pil-image-make-3x3-grid-from-sequence-images
    def image_grid(self, imgs, rows, cols):
        assert len(imgs) == rows * cols

        w, h = imgs[0].size
        grid = Image.new('RGB', size=(cols * w, rows * h))
        grid_w, grid_h = grid.size

        for i, img in enumerate(imgs):
            grid.paste(img, box=(i % cols * w, i // cols * h))
        return grid

    # Produces Circularity-Lz-etc style plots from Sarah Sofie Lovdals thesis.
    # Fine-tuned for our data. Adjust for yours, if you are using this.
    def sofie(self, table):
        # Set up each individua (sub)plot
        axes_xy = [["Lz", "E"], ["Lz", "L"], ["circ", "Lz"], ["L", "E"], ["circ", "L"], ["circ", "E"]]
        lims = {
            "L": [0, 15000],
            "Lz": [-15000,15000],
            "E": [-180000, 100000],
            "circ": [-0.75, 0.75]
        }
        labels = {
            "L": "L",
            "Lz": "Lz",
            "E": "E",
            "circ": r'$\eta$'
        }
        images = []
        for combo in axes_xy:
            xs = table[combo[0]]
            ys = table[combo[1]]
            cmap = cm.get_cmap("hot")
            cmap.set_under('w')
            fig = plt.figure(figsize=(7.5, 5))
            ax = fig.add_subplot(1, 1, 1, projection='scatter_density')

            # Make the norm object to define the image stretch
            #norm = ImageNormalize(vmin=0., vmax=8, stretch=LogStretch())
            axis = ax.scatter_density(xs, ys, cmap=cmap, vmin=.01, vmax=8)#, norm=norm)

            # ax.set(xlim=[-0.1,0.1],
            #       ylim=[-1000,1000])

            fig.colorbar(axis, label='density of stars')

            ax.set(xlim=lims[combo[0]],
                   ylim=lims[combo[1]])
            ax.set(xlabel=labels[combo[0]],
                   ylabel=labels[combo[1]])

            ax.grid(True, which='major', alpha=0, linestyle='dotted')  # Enable grids on subplot
            ax.grid(True, which='minor', alpha=0, linestyle='dotted')
            # ax.set(xlim=lims,
            #       ylim=lims)
            ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                           top=True)
            ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                           top=True)

            image = combo[0] + "_" + combo[1] + ".png"
            images.append(image)
            plt.savefig(image, dpi=300)


        # Clip the images together, now.
        # https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python
        list_im = images
        imgs = [PIL.Image.open(i) for i in list_im]
        grid = self.image_grid(imgs, rows=2, cols=3)
        grid.save(windows_directories.imgdir + "\\" + 'sixplot.jpg')
        for image in list_im:
            try:
                os.remove(image)
            except:
                print("Failed removing image in graphtuils/sofie", image)
                time.sleep(9999999999)

    # Orbit graphing. Take an orbit, a data table, and produce a menagerie of plots demonstrating the fit.
    def orbiplot(self, table, clustering, clust_to_fit, orbit, integration_time, number_of_steps, directory, unique_text, method='dopr54_c'):

        # Get clustering
        truefalse = [True if d == clust_to_fit else False for d in clustering]
        data_to_fit = table[truefalse]

        # Define Orbigist
        import energistics
        orbigist = energistics.orbigistics()

        # Integrate the orbit
        forward = copy.deepcopy(orbit)
        backward = copy.deepcopy(orbit)

        # Do integrals and append these orbit objects to lists
        forward.integrate(t=np.linspace(0, integration_time, number_of_steps)*u.yr,
                          pot=orbigist.pot,
                          method=method)
        backward.integrate(t=-1 * np.linspace(0, integration_time, number_of_steps)*u.yr,
                           pot=orbigist.pot,
                           method=method)

        # Obtain Galactocentric Cylindricals for the model
        orbits = np.concatenate([forward.getOrbit(), backward.getOrbit()], axis=0)
        R, vR, vT, z, vz, phi = orbits.T
        R *= orbigist.rovo[0]
        vR *= orbigist.rovo[1]
        vT *= orbigist.rovo[1]
        z *= orbigist.rovo[0]
        vz *= orbigist.rovo[1]
        model_orbits = [R, vR, vT, z, vz, phi]

        # Also find it for the data
        data_to_fit = orbigist.converter.nowrite_GAL_to_GALCENT(data_to_fit)
        Rdata, vRdata, vTdata, zdata, vzdata, phidata = orbigist.get_leftgalpy(data_to_fit)
        data_orbits = [Rdata, vRdata, vTdata, zdata, vzdata, phidata]

        # Define the rows/columns
        labels = ["R / kpc", "vR / kms" + r'$^{-1}$', "vT / kms" + r'$^{-1}$', "z / kpc",
                  "vz / kms" + r'$^{-1}$', r'$\phi$' + " / rad"]
        saveid = ["R", "vR", "vT", "z", "vz", "phi"]

        # Get permutations for plotting (unique permutations)
        """
        6 unique
        2 in a combo
        6C2 is 6!/(4!*2!) = 15
        Hence a 5x3 graph.
        """
        num_permutations = list(itertools.combinations(np.arange(0,6,1), 2))

        # The plot saveid list for the permutations
        plotids = []

        # Create each individual plot
        for permutation in num_permutations:

            # Set up figure
            fig, ax = plt.subplots(nrows=1,ncols=1)

            # Plot the model
            ax.scatter(model_orbits[permutation[0]], model_orbits[permutation[1]], color='red', s=0.3, marker='x')

            # Plot the data
            ax.scatter(data_orbits[permutation[0]], data_orbits[permutation[1]], color='black', s=0.9, marker='x')

            # Generate the range for the plotting and set lims
            xlength = np.max(data_orbits[permutation[0]]) - np.min(data_orbits[permutation[0]])
            ylength = np.max(data_orbits[permutation[1]]) - np.min(data_orbits[permutation[1]])
            xlim = [np.min(data_orbits[permutation[0]]) - 0.05 * xlength, np.max(data_orbits[permutation[0]]) + 0.05 * xlength]
            ylim = [np.min(data_orbits[permutation[1]]) - 0.05 * ylength, np.max(data_orbits[permutation[1]]) + 0.05 * ylength]
            ax.set(xlim=xlim, ylim=ylim)

            # Generate the legend
            legend_elements = [Patch(edgecolor='gray', facecolor='red', label='Fit'),
                               Patch(edgecolor='gray', facecolor='black', label='Data')]
            plt.legend(handles=legend_elements, loc='upper right')

            # Set the label
            ax.set(xlabel=labels[permutation[0]], ylabel=labels[permutation[1]])
            plt.grid(which='major', color='pink')

            # Generate save ID
            this_save = directory + "\\" + saveid[permutation[0]] + "_" + \
                        saveid[permutation[1]] + "_" + unique_text + ".png"
            plotids.append(this_save)

            # Save the figure
            plt.savefig(this_save, dpi=150)
            plt.close()

        # Clip the images together, now.
        # https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python
        list_im = plotids
        imgs = [PIL.Image.open(i) for i in list_im]
        grid = self.image_grid(imgs, rows=3, cols=5)
        grid.save(directory + "\\" + 'orbifitplot.jpg')

    # Create an energy plot, given a table, and a set of clusters, with colours for each cluster/etc.
    def energy_plot(self, table, clustering, clusters, texts, x_lim, y_lim, savepath=None, mcmillan=False, kpcmyr=False):

        # Energy/etc
        import energistics
        if mcmillan == False:
            energist = energistics.fast_energistics()
            table = energist.default_E_c(table)
        else:
            # Calculate energy via McMillan 2017
            energist = energistics.energistics()
            table = energist.pot_eval(table)


        # Convert units if necessary
        if kpcmyr == True:
            table['E'] = (table['E']*(u.km**2/u.s**2)).to(u.kpc**2/u.myr**2)
            table['Lz'] = (table['Lz']*(u.kpc*u.km/u.s)).to(u.kpc**2/u.myr)

        # Create axis
        fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(10,10))

        # Legend elements
        legend_elements = []

        # Plot all the data
        axs.scatter(table['Lz'], table['E'], marker='o', color='gray', s=8)
        legend_elements.append(Patch(edgecolor='black', facecolor='gray', label='All Data'))
        c = ["red",
             "darkorange",
             "lime",
             "cyan",
             "blue",
             "fuchsia",
             "darkviolet",
             "dodgerblue",
             "chocolate"]
        if len(c) < len(clusters):
            randoms = np.random.random((len(clusters) - len(c), 3))
            for i in randoms:
                c.append(i)

        # For each clustering, plot the energy/etc
        for num,cluster in enumerate(clusters):

            retain = [True if d == cluster else False for d in clustering]
            cluster_table = table[retain]

            cluster_E, cluster_Lz = cluster_table['E'], cluster_table['Lz']

            scat = axs.scatter(cluster_Lz, cluster_E, color=c[num], marker='x', s=16)

            legend_elements.append(Patch(edgecolor='black',
                                         facecolor=c[num],
                                         label=texts[num]))

        # Grid-up
        plt.legend(handles=legend_elements, loc='upper right')
        axs.grid(which='major', color='pink')

        if kpcmyr == False:
            plt.xlabel(r'$\textrm{L}_\textrm{z}\textrm{ / }\textrm{kpc}\cdot\textrm{km}\cdot\textrm{s}^{-1}$')
            plt.ylabel(r'$\textrm{E / km}^2\cdot\textrm{s}^{-2}$')
            axs.set(xlim=x_lim, ylim=y_lim)
        else:
            plt.xlabel(r'$\textrm{L}_\textrm{z}\textrm{ / }\textrm{kpc}^2\cdot\textrm{Myr}^{-1}$')
            plt.ylabel(r'$\textrm{E / kpc}^2\cdot\textrm{Myr}^{-2}$')

            y_lim = (np.array(y_lim)*(u.km**2/u.s**2)).to(u.kpc**2/u.myr**2).value
            x_lim = (np.array(x_lim)*(u.kpc*u.km/u.s)).to(u.kpc**2/u.myr).value

            axs.set(xlim=x_lim, ylim=y_lim)

        if savepath != None:
            try:
                plt.savefig(savepath, dpi=300)
            except:
                plt.savefig(savepath)
        plt.show()