import multiprocessing

import warnings

import numpy as np
from pandas import DataFrame
import ascii_info, windows_directories
from ascii_info import n_carlo
import pickle

import hdfutils
import windows_multiprocessing
from energistics import orbifitter

"""
Monte-carlo orbit fit error estimation. 
Estimate purely based on 200 sets instead of 1,000 due to time constraints (it takes a really fucking long time to run!) 
"""


# Suppress Loky Warnings
warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

# Holder for whether to generate orbits or not. If False, use old orbits.
should_run = True

if __name__ == "__main__":

    # The writer
    writer = hdfutils.hdf5_writer(windows_directories.datadir,
                                  ascii_info.finetune_asciiname)

    # Decide group and the clusters to cluster
    group = ascii_info.fullgroup

    # Take data greattable and obtain clusters to orbifit (length =/= 0)
    clusters_to_orbifit = []
    table = writer.read_table(ascii_info.fullgroup,ascii_info.set_raw)
    prelim_clust = table['prelim_clust']
    set = list(set(prelim_clust))
    for clust in set:

        if clust not in clusters_to_orbifit:

            if clust != -1:

                size_fraction = len(np.where(prelim_clust==clust)[0])/len(prelim_clust)

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if size_fraction < 0.1:

                    clusters_to_orbifit.append(clust)

    # The saveids
    saveids = ascii_info.finetune_orbifit_maindata_saveids

# Run only if name is main.
if __name__ == "__main__" and should_run == True:

    tables = []

    try:

        with open(windows_directories.datadir + "\\finetune_subtableorbits.txt", 'rb') as f:
            tables = pickle.load(file=f)
            tables = tables[0:n_carlo]
            print("Loaded")

    except:

        # Grab the table
        table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)

        # Grab the elements-of-interest
        table = table[['l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost']]

        # Run the pool!
        pool = multiprocessing.Pool(8)
        tables = pool.map(windows_multiprocessing.do_monte_table, [table for i in range(n_carlo)])
        pool.close()

        # Try to save the table
        with open(windows_directories.datadir + "\\finetune_subtableorbits.txt", 'wb') as f:
            pickle.dump(file=f, obj=tables)

    # Parameters for multiprocessing
    parameterss = []

    # Zip the tables up ready for service. Integrate for 1 Gyr
    for tab, saveid in zip(tables, saveids):
        parameterss.append([[group,saveid],
                             tab,
                             clusters_to_orbifit,
                             1000,
                             0.5e9,
                             500])

    # Run the pool!
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.finetune_do_orbifit_maindata, parameterss)
    pool.close()

    # Set the ran
    should_run = False

# Only do this if has_ran is true (get the errors/means/etc on the orbits: orbistatistics.
if should_run == False and __name__ == "__main__": # False:

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
    list_of_stats = pool.map(windows_multiprocessing.finetune_maindata_do_orbistatistics, clusters_to_orbifit)
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

    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.finetune_asciiname)
    writer.write_df("greatfit_monte_orbistatistics_maindata", "data", df)