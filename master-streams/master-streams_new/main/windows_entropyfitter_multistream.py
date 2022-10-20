import copy
import multiprocessing
import pickle
import warnings
import numpy as np
from astropy.table import Table

import ascii_info_new
import hdfutils
import windows_directories_new
import windows_multiprocessing_new
from energistics_new import orbigistics

warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

# Holder for whether to generate orbits or not. If False, use old orbits.
should_run = True
should_load = True

# Decide monte points
nmonte = 100000

# Decide group and the clusters to cluster
group = ascii_info_new.fullgroup
clusts_to_seventropyfit = [[1,3,4,8]]
n_carlo = 50
ranges = np.array([[0.5,2],
                   [0.5,2],
                   [0.5,2]]) # Md Mnfw Cnfw


# Run only if name is main.
if __name__ == "__main__" and should_run == True:

    # Grab the table
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                                  ascii_info_new.fullset)

    # Load in the "mean" percenttable and map
    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    membership_table = writer.read_table(ascii_info_new.fullgroup, "percent_table_greatfitted")
    clustering = membership_table['probable_clust']

    # Grab the elements-of-interest
    table = table[['l','b','dist','dmu_l','dmu_b','vlos','edist','edmu_l','edmu_b','evlost']]

    # Generate n_carlo tables like table
    tables = []

    # Grab orbigist
    orbigist = orbigistics()

    if should_load == True:
        try:
            with open(windows_directories_new.datadir + "\\subtableorbits.txt", 'rb') as f:
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
            with open(windows_directories_new.datadir + "\\subtableorbits.txt", 'wb') as f:
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
        with open(windows_directories_new.datadir + "\\subtableorbits.txt", 'wb') as f:
            pickle.dump(obj=tables, file=f)

            # Parameters for multiprocessing

    parameters = []
    for stable in tables:
        parameters.append([stable, clustering, clusts_to_seventropyfit, ranges, nmonte])

    # Run the pool!
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing_new.do_seventropy, parameters)
    pool.close()

    # Combine all the pointtables for averages/errors  "clust", "M_d", "M_nfw", "c_nfw" returned.
    mds, mnfws, cnfws = np.array(results).T
    md_mean, mnfw_mean, cnfw_mean = np.mean(mds, axis=0), np.mean(mnfws, axis=0), np.mean(cnfws, axis=0)
    mdstd, mnfwstd, cnfwstd = np.std(mds, axis=0), np.std(mnfws, axis=0), np.std(cnfws, axis=0)

    # Save the table for the sake of effort
    dump_array = [mds, mnfws, cnfws]
    with open(windows_directories_new.datadir + "\\multistream_dump_onlyhalo.txt", "wb") as f:
        pickle.dump(obj=dump_array, file=f)

    # Set up table
    save_array = np.array([md_mean, mdstd, mnfw_mean, mnfwstd, cnfw_mean, cnfwstd])

    # Save the table
    savetable = Table(save_array, names=["M_d", "M_d_err", "M_nfw", "M_nfw_err", "c_nfw", "c_nfw_err"])
    writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
    writer.write_table("entropy_fit_witherrors", "multicluster_threeD_onlyhalo", savetable)