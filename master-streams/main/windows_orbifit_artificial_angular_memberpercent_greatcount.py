import copy
import multiprocessing
import os
import time
import warnings
from pathlib import Path

import pandas
from astropy.table import Table
from pandas import DataFrame

import ascii_info, windows_directories
from ascii_info import n_carlo, orbifit_saveids
import pickle
import numpy as np
import graphutils
import hdfutils
import galcentricutils
import windows_multiprocessing
from energistics import orbifitter

"""
Monte-carlo orbit fit error estimation. 
Estimate purely based on 200 sets instead of 1,000 due to time constraints (it takes a really fucking long time to run!) 
"""


# Suppress Loky Warnings
warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

# Holder for whether to generate orbits or not. If False, use old orbits.
should_run = False

# Decide group and the clusters to cluster
group = ascii_info.fullgroup
clusters_to_orbifit = ascii_info.clusters_to_orbifit

# The saveids
saveids = orbifit_saveids

# Run only if name is main.
if __name__ == "__main__" and should_run == True:

    # Grab the table
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)

    # Grab the elements-of-interest
    table = table[['l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost']]

    # Generate n_carlo tables like table
    tables = []

    # Generate error'd tables
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

        # Append
        tables.append(subtable)

    # Parameters for multiprocessing
    parameterss = []

    # Zip the tables up ready for service. Integrate for 1 Gyr
    for tab, saveid in zip(tables, saveids):
        parameterss.append([[group,saveid],
                             tab,
                             clusters_to_orbifit,
                             2000,
                             1e9,
                             2000])

    # Run the pool!
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.do_orbifit, parameterss)
    pool.close()

    # Set the ran
    should_run = False

# Only do this if has_ran is true (get the errors/means/etc on the orbits: orbistatistics.
if should_run == False: # False:

    # The fitter
    orbifit = orbifitter()

    # Set up lists to hold data.
    eesTT, meaneeTT, stdeeTT, \
    periggsTT, meanpgTT, stdpgTT, \
    EEsTT, meanEETT, stdEETT, \
    LzsTT, meanLzTT, stdLzTT, \
    apoggsTT, meanapogTT, stdapogTT = [],[],[],\
                                      [],[],[],\
                                      [],[],[],\
                                      [],[],[],\
                                      [],[],[],

    # For cluster...
    for clust_to_fit in clusters_to_orbifit:

        # Specify savedir/savename and make path, for images
        savedir = windows_directories.imgdir + "\\orbit_fitting_variables" + "\\" + str(clust_to_fit)
        save_unique = str(clust_to_fit)
        try:
            os.mkdir(savedir)
        except:
            pass

        # Get orbits for this cluster
        orbits = []

        # In n_carlo
        for saveid in saveids:
            # Load it
            with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                      saveid + "_fitted_orbit_" + str(clust_to_fit) + ".txt", "rb") as f:
                clust_fitted = pickle.load(file=f)
                orbits.append(clust_fitted)

        # Run stats
        ees, meanee, stdee, \
        periggs, meanpg, stdpg, \
        EEs, meanEE, stdEE, \
        Lzs, meanLz, stdLz, \
        apoggs, meanapog, stdapog = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)

        # Append the stats
        eesTT.append(ees), meaneeTT.append(meanee), stdeeTT.append(stdee), \
        periggsTT.append(periggs), meanpgTT.append(meanpg), stdpgTT.append(stdpg), \
        EEsTT.append(EEs), meanEETT.append(meanEE), stdEETT.append(stdEE), \
        LzsTT.append(Lzs), meanLzTT.append(meanLz), stdLzTT.append(stdLz), \
        apoggsTT.append(apoggs), meanapogTT.append(meanapog), stdapogTT.append(stdapog)

        # Produce various plots
        twod_graph = graphutils.twod_graph()
        twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
        twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
        twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_energy")
        twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_Lz")
        twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Set up column data
    all_data = [clusters_to_orbifit,
                meaneeTT, meanpgTT, meanEETT, meanLzTT, meanapogTT,
                stdeeTT, stdpgTT, stdEETT, stdLzTT, stdapogTT,
                eesTT, periggsTT, EEsTT, LzsTT, apoggsTT]
    colss = ['cluster',
             'e', 'peri', 'E', 'Lz', 'apo',
             'e_std', 'peri_std', 'E_std', 'Lz_std', 'apo_std',
             'e_data', 'peri_data', 'E_data', 'Lz_data', 'apo_data', ]

    # All done- save the df
    df = DataFrame(columns=colss)
    for data, col in zip(all_data, colss):
        df[col] = data


    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    writer.write_df("greatfit_monte_orbistatistics", "data", df)