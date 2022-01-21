import time

import numpy as np
import scipy
from astropy.table import vstack, Table
from numpy import argmin
import pickle
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories

"""
Monte-Carlo for a table [group,set]: for multiprocessed Monte-Carlo simulation. Default 4 workers.
"""
# Default monte-carlo parameters for the multiprocessed table. Change externally if needed.
n_monte = 5
sourcecoord = "solar_info.dat"
def do_monte_table(groupset):
    group,set = groupset
    # Grab table
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(group, set)
    # Set up Monte
    sourcedir = windows_directories.sourcedir
    monte = galcentricutils.monte_angular()
    monte.galdefine(sourcedir, sourcecoord)
    # Run Monte on table
    table = monte.table_monte(table,n_monte)
    # Save the table: watch out for if file is already being written to (retry if it fails.)
    while True:
        try:
            writer.write_table(group, set, table)
            break
        except:
            print("Write Conflict. Sleeping")
            time.sleep(5)
            continue
# do_monte but with covariance matrices: note that we're now dealing with PANDAS DATAFRAMES and NOT ASTROPY TABLES
def do_covmonte_table(groupsetsaveset):
    #print("Covmontecarloing' ", groupsetsaveset)
    group,astropyset,pandaset = groupsetsaveset
    # Grab table
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(group, astropyset)
    # Set up Monte
    sourcedir = windows_directories.sourcedir
    monte = galcentricutils.monte_angular()
    monte.galdefine(sourcedir, sourcecoord)
    # Run Monte on table
    df = monte.table_covmonte(table,n_monte)
    # Save the table: watch out for if file is already being written to (retry if it fails.)
    while True:
        try:
            print("writing ", groupsetsaveset)
            writer.write_df(group, pandaset, df)
            break
        except:
            print("Write Conflict Pandas. Sleeping", groupsetsaveset)
            time.sleep(np.random.randint(0,5))
            continue

    return "done"

    # Done
    #print("Done with ", groupsetsaveset)
# duplimonte for each of our tables. Dumps all pickled results inside subdirectors datadir/duplimonte/group.
# use pickle.load to load the file.
def do_duplimonte_table(groupsetm):
    print(groupsetm)
    group, set, m = groupsetm
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    df = writer.read_df(group, set)
    list_of_columns = galcentricutils.monte_angular().panda_duplimonte(df, m)
    # Generate sets for list_of_tables
    indices = [("L_{}.txt").format(d) for d in [str(d) for d in range(m)]]
    # Pickle and dump all the tables to file
    for column, savedex in zip(list_of_columns, indices):
        with open(windows_directories.duplimontedir + "\\" + group + "\\" + savedex, 'wb') as f:
            pickle.dump(column, f)

"""
Example clustering processing to try and guess good DBS parameters for the circle fitting.
"""
def do_dbstest(paramarray):
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    full_data = writer.read_table(ascii_info.fullgroup, ascii_info.fullset)
    cluster = galcentricutils.cluster()
    cluster.dbs(full_data, *paramarray)

"""
Basic multi-processed KMeans routine for great circle counts, 
- generate grid in theta/phi 
- great circle count each grid
- kmeans each great circle table 
- get a score for how good the kmeans is 
- save html for each great circle plot alongside returning the GCC parameter and GCC score 
- decide a score to use: there's sklearn score default, etc.
"""
# Deprecated
def do_clustering(vararray):
    tableinfo, gcc_params, clustering_params = vararray
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(*tableinfo)
    # Great Circle the table
    table = galcentricutils.greatcount().gcc_table(table, *gcc_params)
    # Cluster the table
    clustered_table = galcentricutils.cluster().gaussmix(table,*clustering_params)
# Deprecated.
def do_gaussaicsbics(vararray):
    tableinfo, gcc_params, clustering_params = vararray
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(*tableinfo)
    table = galcentricutils.greatcount().gcc_table(table, *gcc_params)
    galcentricutils.cluster().gaussaicsbics(table,*clustering_params)

# Various confidence map routines for multiprocessing
"""
Feed dataset information in alongside other info below: apply bonferroni correction, too (see: notes.) 
Take a given pole from the grid, its nphi, and the number of phi regions, and a probability P
Split the GCC generated from this point into the nphi regions
Get the significance of each region using scipy binomial CDF and require it to be less than the bon-corrected lim
Select the lowest significance for this GCC: use this for the final feed-back (maximum confidence = min significance)
Feed back a "True" or "False" alongside the raw calculated CDF for each pole, then save the final result. 
We're using the standard "alpha=0.05" as is typical in statistics for this confidence testing kinda thing
"""
# Note that P is the microscopic probability for each segment of getting a star
# group, set, theta, phi, dtheta, radmin, n_phi, P = vararray
# tableinfo, gccsplit_params, P = vararray
def do_gccsplit_confidence(vararray):
    # Get data
    tableinfo, gccsplit_params, P, savedex = vararray
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(*tableinfo)

    # Split table up into nphi regions
    theta, phi, index, dtheta, radmin, n_phi = gccsplit_params
    splittables, radnumber = galcentricutils.greatcount().gcc_splittable(table, theta, phi, dtheta, radmin, n_phi)

    """
    # Debugging code.
    for num, splittable in enumerate(splittables):
        splittable['k_index'] = [num for d in splittable['x']]
    mastersplit = vstack(splittables)
    graphutils.threed_graph().xmeans_L(mastersplit, savedex, True) """

    sigs = []
    # For each region get the significance
    for splittable in splittables:
        n = len(splittable)
        cdf = scipy.stats.binom.cdf(k=n, n=radnumber, p=P) # below n, not above
        sig = 1 - cdf
        sigs.append(sig)

    # Get the (n) for the phi_region that gives minimum significance/maximum confidence
    minsig = argmin(sigs)
    n_of_splittable_for_minsig = len(splittables[minsig])

    # Pick the minimum significance / maximum confidence (uncorrected)
    minsig = min(sigs)


    # Return maximum significance for this point.thetas, phis, indices, sigmins = results
    return [theta, phi, index, n_of_splittable_for_minsig, minsig]

# Do a quick hdbscan clustering using the provided parameters: designed for duplimonte directory structure.
# Arrayinfo should be [group, saveid] for the duplimonte: minpar as in cluster3d()
def do_hdbscan(arrayinfominpar):
    arrayinfo, minpar = arrayinfominpar
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".txt", 'rb') as f:
        data = pickle.load(f)
    data = np.array(data)
    clustered = galcentricutils.cluster3d().listhdbs(data, minpar)
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".cluster.txt", 'wb') as f:
        pickle.dump(obj=clustered, file=f)
    # Generate a HTML for this clustering under imgdir\\ + "kmeans_html\\duplimonte_kmeanshtml\\" + group_saveid.html
    savedexdir = "kmeans_html\\duplimonte_kmeanshtml\\" + arrayinfo[0] + "\\" + arrayinfo[1]
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, False)



