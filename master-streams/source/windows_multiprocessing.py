import time
import numpy as np
import scipy
from astropy.table import vstack, Table
from matplotlib import pyplot as plt
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
    print("Covmontecarloing' ", groupsetsaveset)
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
    # Do a quick clean on the table to remove dirty L values whose ratio to sigma exceeds n, here n=3.
    # TODO: Use or not? Try without.
    # df = galcentricutils.cluster3d().L_clean(df, 3)
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
# use pickle.load to load the file. Generates artificial angular momenta sets.
def do_duplimonte_table(groupsetm):
    print(groupsetm)
    group, set, m = groupsetm
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    df = writer.read_df(group, set)
    list_of_columns, list_of_columns_2, list_of_columns_3, list_of_columns_4 = galcentricutils.monte_angular().panda_duplimonte(df, m)
    # Generate sets for list_of_tables
    indices, indices_2, indices_3, indices_4 = [("L_{}.txt").format(d) for d in [str(d) for d in range(m)]], \
                                               [("L4D_{}.txt").format(d) for d in [str(d) for d in range(m)]], \
                                               [("LE_{}.txt").format(d) for d in [str(d) for d in range(m)]], \
                                               [("LXYZ_{}.txt").format(d) for d in [str(d) for d in range(m)]]
    # Pickle and dump all the tables to file for LxLyLz
    for column, savedex in zip(list_of_columns, indices):
        with open(windows_directories.duplimontedir + "\\" + group + "\\" + savedex, 'wb') as f:
            pickle.dump(column, f)
    # Pickle and dump all the tables to file for LLzECirc
    for column, savedex in zip(list_of_columns_2, indices_2):
        with open(windows_directories.duplimontedir + "\\" + group + "\\" + savedex, 'wb') as f:
            pickle.dump(column, f)
    # Pickle and dump all the tables to file for LxLyLzE
    for column, savedex in zip(list_of_columns_3, indices_3):
        with open(windows_directories.duplimontedir + "\\" + group + "\\" + savedex, 'wb') as f:
            pickle.dump(column, f)
    # Pickle and dump all the tables to file for LxLyLzXYZ
    for column, savedex in zip(list_of_columns_4, indices_4):
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
# Should be given a list of angular momenta [[lx1,ly1,lz1],[lx2...]...]
def do_hdbscan(arrayinfominpar):
    # Do the clustering
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
# Should be given a list of angular momenta [[L1 Lz1 E1 circ1],[...]...]
def do_hdbscan_L4D(arrayinfominpar):
    arrayinfo, minpar = arrayinfominpar
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".txt", 'rb') as f:
        data = pickle.load(f)
    data = np.array(data)

    # We need to normalize all the variables now, versus angular momentum L (zeroth.)
    """
    Just normalize to the means. See hdbscan_LE for the haphazard way. 
    """
    data = data.T
    data[0] = data[0]/np.mean(np.abs(data[0]))
    data[1] = data[1]/np.mean(np.abs(data[1]))
    data[2] = data[2]/np.mean(np.abs(data[2]))
    data = data.T
    clustered = galcentricutils.cluster3d().listhdbs(data, minpar)
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".cluster.txt", 'wb') as f:
        pickle.dump(obj=clustered, file=f)
    # Generate a HTML for this clustering under imgdir\\ + "kmeans_html\\duplimonte_kmeanshtml\\" + group_saveid.html
    savedexdir = "kmeans_html\\duplimonte_kmeanshtml\\" + arrayinfo[0] + "\\" + arrayinfo[1] + "_WITH_SOFIELOVDAL"
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, False, False)
# Should be given a list of angular momenta [[lx1,ly1,lz1, E1],[lx2...,E2]...]
def do_hdbscan_LE(arrayinfominpar):
    arrayinfo, minpar = arrayinfominpar
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".txt", 'rb') as f:
        data = pickle.load(f)
    data = np.array(data)

    # An important step: we should rescale the energy to match the angular momentum in terms of dimensionality.
    """
    Most energy is 0->-250,000 (10^5)
    Most momenta is -3,500->3500  (10^3) 
    Divide energy by around 100 to make it of equal importance.
    """
    data = data.T
    data[3] *= 1e-2
    data = data.T
    clustered = galcentricutils.cluster3d().listhdbs(data, minpar)
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".cluster.txt", 'wb') as f:
        pickle.dump(obj=clustered, file=f)
    # Generate a HTML for this clustering under imgdir\\ + "kmeans_html\\duplimonte_kmeanshtml\\" + group_saveid.html
    savedexdir = "kmeans_html\\duplimonte_kmeanshtml\\" + arrayinfo[0] + "\\" + arrayinfo[1] + "_WITH_ENERGY"
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, False, False)

# Second round of hdbscan, with greatfit.
def do_hdbscan_greatfit(arrayinfominpar):

    # Set up parameters and clusters_to_cluster
    arrayinfo, minpar, clusters_to_cluster, gccwidths = arrayinfominpar

    # Load in the data
    with open(windows_directories.duplimontedir + "\\" + arrayinfo[0] + "\\" + arrayinfo[1] + ".txt", 'rb') as f:
        data = pickle.load(f)
    data = np.array(data) # as vectors, i.e. 6 columns and len(data) rows. Lx Ly Lz x y z.

    # Set up position data
    data_pos = data[:, 3:6]

    # Angular data
    data_L = data[:, 0:3]

    # Grab greatcircle data
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)

    # If you want the ones that have been "memberpercented" (i.e. limited by surviving monte-carlo)
    #greattable = writer.read_table("greatcircles_preliminary", "greatcircles") # theta phi clust_to_try

    # If you want the general ones (non-monte-carlo-refined ones)
    greattable = writer.read_table("greatcircles_nonmemberpercent_preliminary", "greatcircles") # theta phi clust_to_try


    # Set up greatcount and clust3d and comparitor
    greatcount = galcentricutils.greatcount()
    clust3d = galcentricutils.cluster3d()
    compclust = galcentricutils.compclust()

    # Optionally graph
    grapher = graphutils.threed_graph()
    #specgrapher = graphutils.spec_graph()

    # Load in the average data, for the sake of "remapping" this data (saves us some effort later.)
    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
        mean_cluster = pickle.load(f)

    # Get the current "minimum_excess_index" for unique labelling of the cluster excess
    current_excess_value = np.max(mean_cluster) + 1

    # Identify the indices in "clust_to_try" that match our clusters.
    greatcircle_clusterindices = []
    for cluster in clusters_to_cluster:
        for num, greatcluster in enumerate(greattable['clust_to_try']):
            if greatcluster == cluster:
                greatcircle_clusterindices.append(num)

    # For each cluster in clusters_to_cluster...
    for cluster, great_index, gccwidth in zip(clusters_to_cluster, greatcircle_clusterindices, gccwidths):

        # Grab the great circle
        theta, phi = greattable[great_index]['theta'], greattable[great_index]['phi']

        # Cut the data by the greatcircle (angular momentum-ly)
        truefalse, trues = greatcount.gcc_array_retain(data_pos, theta, phi, gccwidth)

        cluster_L = data_L[truefalse]

        # Cluster the greatcircle cut
        clustered_L = clust3d.listhdbs(cluster_L, minpar)

        # Need to bring clustered_L back to the original length/shape (we only clustered a subsection of meandata)
        # Replace all greatfit-removed points by -1 (noise)
        clustered = np.ones(len(truefalse), dtype=int)
        clustered *= -1
        clustered[trues] = clustered_L

        # Match to the mean clustering
        remap = compclust.compclust_multilabel(mean_cluster, clustered, current_excess_value)

        # Save the remap
        with open(windows_directories.duplimontedir + "\\" + ascii_info.fullgroup +
                  "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}.txt").format(cluster), 'wb') as f:
            pickle.dump(obj=remap, file=f)

        # Generate graphs, too. For debug, mainly.
        grapher.kmeans_L_array(cluster_L, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_ang").format(cluster), False, True)
        # Also positions for graphing
        #cluster_pos = np.array(data_pos[truefalse])
        #x,y,z = cluster_pos.T
        #plt.scatter(y,x)
        #plt.show()
        #grapher.xmeans_L_array(cluster_pos, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_pos").format(cluster), False, True)
        #time.sleep(500)
        #specgrapher.array_thetaphi_debug(cluster_pos, theta, phi, remap[truefalse], cluster)
    #time.sleep(5000)

# Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
def do_orbifit(parameters):
    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit
        clust_fitted = fitter.galpy_final_fitting(table, clust_to_fit, iterations, time_to_integrate,
                                                  number_of_steps, try_load=True, graph=False,
                                                  load_fit=False, try_save=False, debug_graph=None)

        # Save it
        with open(windows_directories.orbitsfitdir + "\\" + arrayinfo[0] + "_" +
                  arrayinfo[1] + "_fitted_orbit_" + str(clust_to_fit) + ".txt", "wb") as f:
            pickle.dump(obj=clust_fitted, file=f)

# Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
# Specific for maindata (no membership table applicable- just the maindata clustering.)
def do_orbifit_maindata(parameters):
    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # Grab the clustering
    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
        clustered = pickle.load(file=f)

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit
        clust_fitted = fitter.galpy_fitting_nomemb(table, clustered, clust_to_fit, iterations, time_to_integrate,
                                                   number_of_steps, try_load=True, graph=False,
                                                   load_fit=False, try_save=False, extra_text="orbifit_maindata_run")

        # Save it
        with open(windows_directories.orbitsfitdir + "\\" + arrayinfo[0] + "_" +
                  arrayinfo[1] + "_fitted_orbit_maindata_" + str(clust_to_fit) + ".txt", "wb") as f:
            pickle.dump(obj=clust_fitted, file=f)