import copy
import multiprocessing
import os
import time
import warnings
from pathlib import Path
import pandas
from astropy.table import Table
from pandas import DataFrame
import ascii_info_new, windows_directories_new
from ascii_info_new import n_carlo, orbifit_saveids
import pickle
import numpy as np
import graphutils_new
import hdfutils
import galcentricutils_new
import windows_multiprocessing_new
from energistics_new import orbifitter

"""
Monte-carlo orbit fit error estimation. 
Estimate purely based on 200 sets instead of 1,000 due to time constraints (it takes a really fucking long time to run!) 
"""


# Suppress Loky Warnings
warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

# Holder for whether to generate orbits or not. If False, use old orbits.
should_run = True

# Decide group and the clusters to cluster
group = ascii_info_new.fullgroup
clusters_to_orbifit = ascii_info_new.clusters_to_maindata_orbifit

# The saveids
saveids = ascii_info_new.orbifit_maindata_saveids

# Run only if name is main.
if __name__ == "__main__" and should_run == True:

    tables = []

    try:
        with open(windows_directories_new.datadir + "\\subtableorbits.txt", 'rb') as f:
            tables = pickle.load(file=f)
            tables = tables[0:n_carlo]
            print("Loaded")
    except:

        # Grab the table
        table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                     ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                                      ascii_info_new.fullset)

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
                             1000])

    # Run the pool!
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing_new.do_orbifit_maindata, parameterss)
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
                                      [],[],[]

    # Run the pool!
    pool = multiprocessing.Pool(8)
    list_of_stats = pool.map(windows_multiprocessing_new.do_orbistatistics, clusters_to_orbifit)
    pool.close()

    # Iterate over the list of stats
    for stats in list_of_stats:

        # Section stats
        ees, meanee, stdee, \
        periggs, meanpg, stdpg, \
        EEs, meanEE, stdEE, \
        Lzs, meanLz, stdLz, \
        apoggs, meanapog, stdapog = stats

        # Append the stats
        eesTT.append(ees), meaneeTT.append(meanee), stdeeTT.append(stdee), \
        periggsTT.append(periggs), meanpgTT.append(meanpg), stdpgTT.append(stdpg), \
        EEsTT.append(EEs), meanEETT.append(meanEE), stdEETT.append(stdEE), \
        LzsTT.append(Lzs), meanLzTT.append(meanLz), stdLzTT.append(stdLz), \
        apoggsTT.append(apoggs), meanapogTT.append(meanapog), stdapogTT.append(stdapog)

    # Set up column data
    all_data = [clusters_to_orbifit,
                meaneeTT, meanpgTT, meanEETT, meanLzTT, meanapogTT,
                stdeeTT, stdpgTT, stdEETT, stdLzTT, stdapogTT,
                eesTT, periggsTT, EEsTT, LzsTT, apoggsTT]
    colss = ['cluster',
             'e', 'peri', 'E', 'Lz', 'apo',
             'e_std', 'peri_std', 'E_std', 'Lz_std', 'apo_std',
             'e_data', 'peri_data', 'E_data', 'Lz_data', 'apo_data']

    # All done- save the df
    df = DataFrame(columns=colss)
    for data, col in zip(all_data, colss):
        df[col] = data


    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    writer.write_df("greatfit_monte_orbistatistics", "data", df)