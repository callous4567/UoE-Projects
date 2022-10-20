import copy
import multiprocessing
import pickle
import warnings
import numpy as np
from astropy.table import Table

import ascii_info
import hdfutils
import windows_directories
import windows_multiprocessing
from energistics import orbigistics



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
import ascii_info
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

warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

# Holder for whether to generate orbits or not. If False, use old orbits.
should_run = True
should_load = True

# Decide monte points
nmonte = 100

# Decide group and the clusters to cluster
group = ascii_info.fullgroup
clusters_to_singleentropyfit = [1,3,4,5,6,8,11]
n_carlo = 50
ranges = np.array([[1,1],
                   [0.5,2],
                   [1,1]]) # Md Mnfw Cnfw


# Run only if name is main.
if __name__ == "__main__" and should_run == True:

    # Grab the table
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)

    # Load in the "mean" percenttable and map
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    membership_table = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")
    clustering = membership_table['probable_clust']

    # Grab the elements-of-interest
    table = table[['l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost']]

    # Generate n_carlo tables like table
    tables = []

    # Grab orbigist
    orbigist = orbigistics()

    if should_load == True:
        try:
            with open(windows_directories.datadir + "\\subtableorbits.txt", 'rb') as f:
                tables = pickle.load(file=f)
                tables = tables[0:n_carlo]
                print("Loaded")
        except:
            # Generate error'd tables
            print("Load failed")
            for i in range(n_carlo):

                # Copy the table
                subtable = copy.deepcopy(table)

                # Now, for the rows.
                for row in subtable:

                    # Generate artificial data
                    dist = np.random.default_rng().normal(row['dist'], row['edist'], 1)
                    dist = np.abs(dist) # negatives happen- bad for calculating.
                    dmul = np.random.default_rng().normal(row['dmu_l'], row['edmu_l'], 1)
                    dmub = np.random.default_rng().normal(row['dmu_b'], row['edmu_b'], 1)
                    vlos = np.random.default_rng().normal(row['vlos'], row['evlost'], 1)
                    row['dist'], row['dmu_l'], row['dmu_b'], row['vlos'] = dist[0], dmul[0], dmub[0], vlos[0]

                # Convert galactic coordinates to galactocentric coordinates
                subtable = orbigist.converter.nowrite_GAL_to_GALCENT(subtable)

                # Append
                tables.append(subtable)
            with open(windows_directories.datadir + "\\subtableorbits.txt", 'wb') as f:
                pickle.dump(obj=tables, file=f)

    else:
        # Generate error'd tables
        for i in range(n_carlo):

            # Copy the table
            subtable = copy.deepcopy(table)

            # Now, for the rows.
            for row in subtable:
                # Generate artificial data
                dist = np.random.default_rng().normal(row['dist'], row['edist'], 1)
                dist = np.abs(dist)  # negatives happen- bad for calculating.
                dmul = np.random.default_rng().normal(row['dmu_l'], row['edmu_l'], 1)
                dmub = np.random.default_rng().normal(row['dmu_b'], row['edmu_b'], 1)
                vlos = np.random.default_rng().normal(row['vlos'], row['evlost'], 1)
                row['dist'], row['dmu_l'], row['dmu_b'], row['vlos'] = dist[0], dmul[0], dmub[0], vlos[0]

            # Convert galactic coordinates to galactocentric coordinates
            subtable = orbigist.converter.nowrite_GAL_to_GALCENT(subtable)

            # Append
            tables.append(subtable)
        with open(windows_directories.datadir + "\\subtableorbits.txt", 'wb') as f:
            pickle.dump(obj=tables, file=f)

            # Parameters for multiprocessing

    parameters = []
    for stable in tables:
        parameters.append([stable, clustering, clusters_to_singleentropyfit, ranges, nmonte])

    # Run the pool!
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.do_minentropy, parameters)
    pool.close()

    # Combine all the pointtables for averages/errors  "clust", "M_d", "M_nfw", "c_nfw" returned.
    mds, mnfws, cnfws = [[] for d in clusters_to_singleentropyfit], \
                        [[] for d in clusters_to_singleentropyfit], \
                        [[] for d in clusters_to_singleentropyfit]
    for pointtable in results:
        for i in range(len(clusters_to_singleentropyfit)):
            mds[i].append(pointtable['M_d'][i])
            mnfws[i].append(pointtable['M_nfw'][i])
            cnfws[i].append(pointtable['c_nfw'][i])

    # Get means
    mds, mnfws, cnfws = np.array(mds), np.array(mnfws), np.array(cnfws)
    md_mean, mnfw_mean, cnfw_mean = np.mean(mds, axis=1), np.mean(mnfws, axis=1), np.mean(cnfws, axis=1)
    mdstd, mnfwstd, cnfwstd = np.std(mds, axis=1), np.std(mnfws, axis=1), np.std(cnfws, axis=1)

    # Set up table
    save_array = np.empty(shape=(len(clusters_to_singleentropyfit), 7))
    save_array[:,0] = clusters_to_singleentropyfit
    save_array[:,1] = md_mean
    save_array[:,2] = mdstd
    save_array[:,3] = mnfw_mean
    save_array[:,4] = mnfwstd
    save_array[:,5] = cnfw_mean
    save_array[:,6] = cnfwstd

    # Save the table
    savetable = Table(save_array, names=["cluster", "M_d", "M_d_err", "M_nfw", "M_nfw_err", "c_nfw", "c_nfw_err"])
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    writer.write_table("entropy_fit_witherrors", "per_cluster_haloonly", savetable)