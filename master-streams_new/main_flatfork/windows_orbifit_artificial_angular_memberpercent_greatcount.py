import multiprocessing
import warnings
from pandas import DataFrame
import ascii_info_new, windows_directories_new
from ascii_info_new import n_carlo, orbifit_saveids
import pickle
import hdfutils
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

# Take data greattable and obtain clusters to orbifit (length =/= 0)
if __name__ == "__main__":

    clusters_to_orbifit = []
    total_percent_table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                  "total_percent_table_greatfitted")

    for num,clust in enumerate(total_percent_table['cluster']):

        if clust not in clusters_to_orbifit:

            if clust != -1:

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if total_percent_table[num]['membership_fraction'] < 0.1:

                    clusters_to_orbifit.append(clust)

    # Print it
    print("Orbifitting for ", clusters_to_orbifit)
    import time
    time.sleep(10)

#clusters_to_orbifit = ascii_info_new.flatfork_clusters_to_maindata_orbifit

# The saveids
saveids = orbifit_saveids

# Run only if name is main.
if __name__ == "__main__" and should_run == True:

    # Grab the table
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                  ascii_info_new.fullset)

    # Grab the elements-of-interest
    table = table[['l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost']]

    # Generate n_carlo tables like table
    tables = []

    try:

        with open(windows_directories_new.datadir + "\\flatfork_subtableorbits.txt", 'rb') as f:
            tables = pickle.load(file=f)
            tables = tables[0:n_carlo]
            print("Loaded")

    except:

        # Grab the table
        table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                     ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                      ascii_info_new.fullset)

        # Grab the elements-of-interest
        table = table[['l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost']]

        # Run the pool for the preliminaries
        pool = multiprocessing.Pool(8)
        tables = pool.map(windows_multiprocessing_new.do_monte_table, [table for i in range(n_carlo)])
        pool.close()

        # Try to save the table
        with open(windows_directories_new.datadir + "\\flatfork_subtableorbits.txt", 'wb') as f:
            pickle.dump(file=f, obj=tables)

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

    # Run the pool! For preliminaries first
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing_new.flatfork_do_orbifit, parameterss)
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
    list_of_stats = pool.map(windows_multiprocessing_new.flatfork_do_orbistatistics, clusters_to_orbifit)
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


    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.flatfork_asciiname)
    writer.write_df("greatfit_monte_orbistatistics", "data", df)