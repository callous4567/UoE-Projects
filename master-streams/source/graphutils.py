import copy
import os
import time

import astropy
import numpy as np
from astropy.table import Table, unique, vstack

# Our own stuff.

import galcentricutils
import hdfutils
import windows_directories

# Matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import matplotlib.colors as mcolors
from matplotlib import cm


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

# Holds all graphing utilities for project. Notes:
"""
Units should be as standard.
xyz in kpc
vxvyvz in kms^-1
lxlylz in kpc kms^-1
angles in degrees
"""

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

    # Plot nclust vs. n for a monte-carloed clustering set.
    # Also subplot the PDF for the clustering.
    def nclust_n(self, data, group):
        # The n_clust/clustering plot.
        x = np.arange(0, len(data), 1)
        y = data
        fig, axs = plt.subplots(nrows=1,ncols=2,sharey='row', gridspec_kw={'width_ratios': [3,1]})
        plt.subplots_adjust(wspace=0.1, hspace=0)
        axs[0].scatter(x, y, marker='x',color='black',s=10)
        axs[0].set(xlabel="clustering/index",
               ylabel="n_clust")

        # Generate a PDF, too.
        axs[1].hist(y, bins=50, density=True, orientation='horizontal', color='black')
        plt.suptitle("For " + group)
        plt.show()


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
    def kmeans_L_array(self, array, clusterdata, savedexdir, browser):
        x,y,z = array.T
        fig = px.scatter_3d(x=x, y=y, z=z, color=clusterdata)
        if savedexdir != False :
            save_loc = windows_directories.imgdir + "\\" + savedexdir + ".html"
            fig.write_html(save_loc)
        if browser == True:
            fig.show()

    # Given two clusterings [data1,data2] and label sets [labels1,labels2] plot both
    # Shifts the 2nd labelling by 16 for colours.
    def kmeans_L_array_pair(self, arrays, clusterdatas, savedexdir, browser, outliers=False):
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
    """
    def clust_radec(self, table, clustering, cluster_id=False, vasiliev=False):
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
                tab_x, tab_y = sim_table['ra'], sim_table['dec']
                axs.scatter(tab_x, tab_y, color="gray", s=0.01, marker=".")

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
                tab_x, tab_y = sorted_table['ra'], sorted_table['dec']
                axs.scatter(tab_x, tab_y, s=scale(len(sorted_table)), marker="s")

            # Misc axis things
            axs.set_facecolor("k")
            plt.gca().invert_xaxis()
            axs.grid(True, which='major', alpha=1, linewidth=0.25)
            axs.grid(True, which='minor', alpha=1, linewidth=0.25)
            axs.grid(color="white")
            axs.set(xlabel=r'$\alpha$', ylabel=r'$\delta$')

            # Show
            plt.show()

        else:
            # Grab the ra/dec list for the cluster_id
            table['cluster'] = clustering
            table_by_clust = table.group_by("cluster")
            clustmask = table_by_clust.groups.keys["cluster"] == cluster_id
            masked_table = table_by_clust.groups[clustmask]

            # ICRS-up the thing and do the plot...
            fig, axs = plt.subplots(nrows=1, ncols=1, dpi=300)
            # first do Vasilievs, though.
            if vasiliev==True:
                simdir = windows_directories.datadir + "\\vasiliev\\Sgr_snapshot"
                simfile = "simdata.hdf5"
                writer = hdfutils.hdf5_writer(simdir, simfile)
                sim_table = writer.read_table("Sgr_snapshot", "astrotable")
                tab_x, tab_y = sim_table['ra'], sim_table['dec']
                axs.scatter(tab_x, tab_y, color="gray", s=0.01, marker=".")
            # Then our clustering...
            masked_table = galcentricutils.galconversion().nowrite_GAL_to_ICRS(masked_table)
            tab_x, tab_y = masked_table['ra'], masked_table['dec']
            axs.grid(True, which='major', alpha=1, linewidth=0.25, color="black")  # Enable grids on subplot
            axs.grid(True, which='minor', alpha=1, linewidth=0.25, color="black")
            axs.scatter(tab_x, tab_y, color="red", s=5, marker="s")
            axs.set_facecolor("k")
            plt.gca().invert_xaxis()
            axs.set(xlabel=r'$\alpha$', ylabel=r'$\delta$')
            axs.grid(color="white")
            # We can also add a great circle...
            #gcc = galcentricutils.greatcount().gcc_gen(10000, 45, 5)  # thetas, phis, degrees.
            #axs.scatter(gcc[1], 90 - gcc[0], color="blue", s=0.1)
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
        grid.save('Trifecta.jpg')

