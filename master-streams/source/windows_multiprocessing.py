import os
import time
import numpy as np
from astropy.table import Table
import pickle
import ascii_info
import energistics
import galcentricutils
import graphutils
import hdfutils
import windows_directories
from energistics_constants import M_d, M_b, M_nfw, a_b, a_d, b_d, a_nfw, c_nfw

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
    table['dist'] = np.abs(table['dist'])
    table['edist'] = np.abs(table['edist'])
    table['edmu_l'], table['edmu_b'] = np.abs(table['edmu_l']), np.abs(table['edmu_b'])
    table['evlost'] = np.abs(table['evlost'])
    # Set up Monte
    sourcedir = windows_directories.sourcedir
    monte = galcentricutils.monte_angular()
    monte.galdefine(sourcedir, sourcecoord)
    # Run Monte on table
    df = monte.table_covmonte(table,n_monte)

    # Save the table: watch out for if file is already being written to (retry if it fails.)
    i = 0
    while i < 10:
        try:
            print("writing ", groupsetsaveset)
            writer.write_df(group, pandaset, df)
            break
        except Exception as e:
            print(e)
            print("Write Conflict Pandas. Sleeping", groupsetsaveset)
            time.sleep(np.random.randint(0,5))
            i += 1
            continue

    return "done"

# duplimonte for each of our tables. Dumps all pickled results inside subdirectors datadir/duplimonte/group.
# use pickle.load to load the file. Generates artificial angular momenta sets.
def do_duplimonte_table(groupsetm):
    print(groupsetm)
    group, set, m = groupsetm
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    df = writer.read_df(group, set)
    list_of_columns, list_of_columns_2, list_of_columns_3, list_of_columns_4 = galcentricutils.monte_angular().panda_duplimonte(df, m) # , list_of_columns_2, list_of_columns_3, list_of_columns_4
    # Generate sets for list_of_tables
    indices , indices_2, indices_3, indices_4 = [("L_{}.txt").format(d) for d in [str(d) for d in range(m)]], \
                                               [("L4D_{}.txt").format(d) for d in [str(d) for d in range(m)]], \
                                               [("LE_{}.txt").format(d) for d in [str(d) for d in range(m)]], \
                                               [("LXYZ_{}.txt").format(d) for d in [str(d) for d in range(m)]]
    # Pickle and dump all the tables to file for LxLyLz
    for column, savedex in zip(list_of_columns, indices):
        with open(os.path.join(os.path.join(windows_directories.duplimontedir,group),savedex), 'wb') as f:
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

def do_hdbscan(arrayinfominpar):
    """

    Do a quick hdbscan clustering using the provided parameters: designed for duplimonte directory structure.
    Arrayinfo should be [group, saveid] for the duplimonte: minpar as in cluster3d()
    Should be given a list of angular momenta [[lx1,ly1,lz1],[lx2...]...]

    :param arrayinfominpar: arrayinfo, minpar = arrayinfominpar
    :return:
    """

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

def flatfork_do_hdbscan(arrayinfominpar):

    """

    ### Flatfork Edition! Load in the flatfork clustering and use that. ###

    Do a quick hdbscan clustering using the provided parameters: designed for duplimonte directory structure.
    Arrayinfo should be [group, saveid] for the duplimonte: minpar as in cluster3d()
    Should be given a list of angular momenta [[lx1,ly1,lz1],[lx2...]...]

    :param arrayinfominpar: arrayinfo, minpar = arrayinfominpar
    :return:
    """

    # Grab the flat clusterer
    with open(windows_directories.clusterer_flat_path, "rb") as f:
        clusterer = pickle.load(file=f)

    # Get parameters/etc + load data
    arrayinfo, minpar = arrayinfominpar
    with open(os.path.join(os.path.join(windows_directories.duplimontedir,
                                        arrayinfo[0]),
                           arrayinfo[1] + ".txt"),
              'rb') as f:
        data = pickle.load(f)
    data = np.array(data)

    # Fit-predict the data
    from hdbscan import prediction
    clustered = prediction.approximate_predict(clusterer, data)[0]
    with open(os.path.join(os.path.join(windows_directories.duplimontedir,
                                        arrayinfo[0]),
                           arrayinfo[1] + "_flatfork_" + ".cluster.txt"),
              'wb') as f:
        pickle.dump(obj=clustered, file=f)

    # Generate a HTML for this clustering under imgdir\\ + "kmeans_html\\duplimonte_kmeanshtml\\" + group_saveid.html
    savedexdir = "flatfork_" + "kmeans_html\\duplimonte_kmeanshtml\\" + arrayinfo[0] + "\\" + arrayinfo[1]
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, False)

def finetune_do_hdbscan(arrayinfominpar):

    """

    ### finetune Edition! Load in the finetune clustering and use that. ###

    Do a quick hdbscan clustering using the provided parameters: designed for duplimonte directory structure.
    Arrayinfo should be [group, saveid] for the duplimonte: minpar as in cluster3d()
    Should be given a list of angular momenta [[lx1,ly1,lz1],[lx2...]...]

    :param arrayinfominpar: arrayinfo, minpar = arrayinfominpar
    :return:
    """

    # Get parameters/etc + load data
    arrayinfo, minpar = arrayinfominpar
    with open(os.path.join(os.path.join(windows_directories.duplimontedir,
                                        arrayinfo[0]),
                           arrayinfo[1] + ".txt"),
              'rb') as f:
        data = pickle.load(f)
    data = np.array(data)

    # Set up the flat hdbscan run. 25 works well (See Flatfork- 25/26.) This is more predictable to produce a mean.
    from hdbscan import flat
    clustered = flat.HDBSCAN_flat(X=data,
                                  n_clusters=ascii_info.fine_nclust,
                                  min_cluster_size=minpar[0],
                                  min_samples=minpar[1],
                                  metric='l2',
                                  algorithm='best',
                                  prediction_data=False).labels_

    with open(os.path.join(os.path.join(windows_directories.duplimontedir,
                                        arrayinfo[0]),
                           arrayinfo[1] + "_finetune" + ".cluster.txt"),
              'wb') as f:
        pickle.dump(obj=clustered, file=f)

    # Generate a HTML for this clustering under imgdir\\ + "kmeans_html\\duplimonte_kmeanshtml\\" + group_saveid.html
    savedexdir = "finetune_" + "kmeans_html\\duplimonte_kmeanshtml\\" + arrayinfo[0] + "\\" + arrayinfo[1]
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
    #grapher = graphutils.threed_graph()
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
        # We justify this: only bothering with the stars in the GCC for memberpercenttable, so non-GCC ones can sod off.
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
        #grapher.kmeans_L_array(cluster_L, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_ang").format(cluster), False, True)
        # Also positions for graphing
        #cluster_pos = np.array(data_pos[truefalse])
        #x,y,z = cluster_pos.T
        #plt.scatter(y,x)
        #plt.show()
        #grapher.xmeans_L_array(cluster_pos, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_pos").format(cluster), False, True)
        #time.sleep(500)
        #specgrapher.array_thetaphi_debug(cluster_pos, theta, phi, remap[truefalse], cluster)
    #time.sleep(5000)

# Flatfork/updated rendition
massive_factor = 1.1
def do_updated_hdbscan_greatfit_flatfork(arrayinfominpar):

    # Split input
    group_saveid, best_minpars, clustrange, greatpars = arrayinfominpar

    # Load in the data
    with open(windows_directories.duplimontedir + "\\" + group_saveid[0] + "\\" + group_saveid[1] + ".txt", 'rb') as f:
        data = pickle.load(f)
    data = np.array(data) # as vectors, i.e. 6 columns and len(data) rows. Lx Ly Lz x y z.

    # Set up position data
    data_pos = data[:, 3:6]

    # Angular data
    data_L = data[:, 0:3]

    # Set up greatcount and clust3d and comparitor
    greatcount = galcentricutils.greatcount()
    clust3d = galcentricutils.cluster3d()
    compclust = galcentricutils.compclust()

    # Load in the average data, for the sake of "remapping" this data (saves us some effort later.)
    with open(windows_directories.clusterdir + "\\" + "flatfork_fullgroup_cluster.txt", 'rb') as f:
        clustered = pickle.load(f)

    # Get the current "minimum_excess_index" for unique labelling of the cluster excess
    current_excess_value = np.max(clustered) + 1

    # For each cluster in clusters_to_cluster...
    for cluster, best_minpar, greatpar in zip(clustrange, best_minpars, greatpars):

        # Grab the great circle
        width, theta, phi = greatpar

        # Cut the data by the greatcircle (angular momentum-ly + also from the mean cluster)
        truefalse, trues = greatcount.gcc_array_retain(data_pos, theta, phi, width)
        cut_clustered = clustered[truefalse]
        gcc_data_L = data_L[truefalse]

        # Assuage how much GSE (the largest cluster other than noise) is in this set of data
        is_gse = np.array([True if d == ascii_info.flatfork_GSE_ID else False for d in cut_clustered])
        is_gse = np.sum(is_gse) * (massive_factor)

        # Cluster the greatcircle cut
        trial_clustered = clust3d.listhdbs_updated(gcc_data_L, best_minpar, max_cluster_size=is_gse)

        # Match to the mean clustering
        remap = np.array(compclust.compclust_multilabel_julia(cut_clustered, trial_clustered, current_excess_value),
                         int)

        # Set those within (the Trues of the GCC) to the remap, and set those outside to noise
        trial_clustered = np.ones(len(truefalse), dtype=int)
        trial_clustered *= -1
        trial_clustered[trues] = remap

        # Save the remap
        with open(windows_directories.duplimontedir + "\\" + group_saveid[0] +
                  "\\" + group_saveid[1] + (".flatfork_remap-cluster_{0:.0f}.txt").format(cluster), 'wb') as f:
            pickle.dump(obj=trial_clustered, file=f)

        # Generate graphs, too. For debug, mainly.
        #grapher.kmeans_L_array(cluster_L, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_ang").format(cluster), False, True)
        # Also positions for graphing
        #cluster_pos = np.array(data_pos[truefalse])
        #x,y,z = cluster_pos.T
        #plt.scatter(y,x)
        #plt.show()
        #grapher.xmeans_L_array(cluster_pos, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_pos").format(cluster), False, True)
        #time.sleep(500)
        #specgrapher.array_thetaphi_debug(cluster_pos, theta, phi, remap[truefalse], cluster)
    #time.sleep(5000)
# Finetune
def do_updated_hdbscan_greatfit_finetune(arrayinfominpar):

    # Split input
    group_saveid, best_minpars, clustrange, greatpars = arrayinfominpar

    # Load in the data
    with open(windows_directories.duplimontedir + "\\" + group_saveid[0] + "\\" + group_saveid[1] + ".txt", 'rb') as f:
        data = pickle.load(f)
    data = np.array(data) # as vectors, i.e. 6 columns and len(data) rows. Lx Ly Lz x y z.

    # Set up position data
    data_pos = data[:, 3:6]

    # Angular data
    data_L = data[:, 0:3]

    # Set up greatcount and clust3d and comparitor
    greatcount = galcentricutils.greatcount()
    clust3d = galcentricutils.cluster3d()
    compclust = galcentricutils.compclust()

    # Load in the average data, for the sake of "remapping" this data (saves us some effort later.)
    with open(windows_directories.clusterdir_fine + "\\" + "finetune_fullgroup.cluster.txt", 'rb') as f:
        clustered = pickle.load(f)

    # Get the current "minimum_excess_index" for unique labelling of the cluster excess
    current_excess_value = np.max(clustered) + 1

    # For each cluster in clusters_to_cluster...
    for cluster, best_minpar, greatpar in zip(clustrange, best_minpars, greatpars):

        # Grab the great circle
        width, theta, phi = greatpar

        # Cut the data by the greatcircle (angular momentum-ly + also from the mean cluster)
        truefalse, trues = greatcount.gcc_array_retain(data_pos, theta, phi, width)
        cut_clustered = clustered[truefalse]
        gcc_data_L = data_L[truefalse]

        # Assuage how much GSE (the largest cluster other than noise) is in this set of data
        is_gse = np.array([True if d == ascii_info.finetune_GSE_ID else False for d in cut_clustered])
        is_gse = np.sum(is_gse) * (massive_factor)

        # Cluster the greatcircle cut
        trial_clustered = clust3d.listhdbs_updated(gcc_data_L, best_minpar, max_cluster_size=is_gse)

        # Match to the mean clustering
        remap = np.array(compclust.compclust_multilabel_julia(cut_clustered, trial_clustered, current_excess_value),
                         int)

        # Set those within (the Trues of the GCC) to the remap, and set those outside to noise
        trial_clustered = np.ones(len(truefalse), dtype=int)
        trial_clustered *= -1
        trial_clustered[trues] = remap

        # Save the remap
        with open(windows_directories.duplimontedir + "\\" + group_saveid[0] +
                  "\\" + group_saveid[1] + (".finetune_remap-cluster_{0:.0f}.txt").format(cluster), 'wb') as f:
            pickle.dump(obj=trial_clustered, file=f)

        # Generate graphs, too. For debug, mainly.
        #grapher.kmeans_L_array(cluster_L, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_ang").format(cluster), False, True)
        # Also positions for graphing
        #cluster_pos = np.array(data_pos[truefalse])
        #x,y,z = cluster_pos.T
        #plt.scatter(y,x)
        #plt.show()
        #grapher.xmeans_L_array(cluster_pos, remap[truefalse],"\\" + "cluster_greatcount_debug" + "\\" + arrayinfo[1] + (".remap-cluster_{0:.0f}_pos").format(cluster), False, True)
        #time.sleep(500)
        #specgrapher.array_thetaphi_debug(cluster_pos, theta, phi, remap[truefalse], cluster)
    #time.sleep(5000)

# Monte a table
def do_monte_table(table):

    # Set up rng
    rng = np.random.default_rng()

    # Monte the entire thing (once.)
    dist = rng.normal(table['dist'], table['edist'], len(table['dist']))
    table['dist'] = np.abs(dist)  # negatives happen- bad for calculating.
    table['dmu_l'] = rng.normal(table['dmu_l'], table['edmu_l'], len(table['dmu_l']))
    table['dmu_b'] = rng.normal(table['dmu_b'], table['edmu_b'], len(table['dmu_b']))
    table['vlos'] = rng.normal(table['vlos'], table['evlost'], len(table['vlos']))

    # Return
    return table

# Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
# Uses membership table
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

# Do the preliminary fit for finetune_do_orbifit
def finetune_do_preliminary(parameters):

    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit
        clust_fitted = fitter.finetune_galpy_final_preliminary(table, clust_to_fit, iterations, time_to_integrate,
                                                  number_of_steps, try_load=True, graph=False,
                                                  load_fit=False, try_save=False, debug_graph=None)

# Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
# Uses membership table
def finetune_do_orbifit(parameters):

    """
    Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
        - Uses membership table

    :param parameters: arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters
    :return: pickle-dump
    """

    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit
        clust_fitted = fitter.finetune_galpy_final_fitting(table, clust_to_fit, iterations, time_to_integrate,
                                                  number_of_steps, try_load=True, graph=False,
                                                  load_fit=False, try_save=False, debug_graph=None)

        # Save it
        with open(windows_directories.orbitsfitdir + "\\" + arrayinfo[0] + "_" +
                  arrayinfo[1] + "_finetune_fitted_orbit_" + str(clust_to_fit) + ".txt", "wb") as f:
            pickle.dump(obj=clust_fitted, file=f)


# Do the preliminary fit for flatfork_do_orbifit
def flatfork_do_preliminary(parameters):

    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit
        clust_fitted = fitter.flatfork_galpy_final_preliminary(table, clust_to_fit, iterations, time_to_integrate,
                                                  number_of_steps, try_load=True, graph=False,
                                                  load_fit=False, try_save=False, debug_graph=None)

# Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
# Uses membership table
def flatfork_do_orbifit(parameters):

    """
    Fit an orbit with the energistics orbifitter using Galpy, then save the fit, for the Monte-Carlo Data Iterations
        - Uses membership table

    :param parameters: arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters
    :return: pickle-dump
    """

    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit
        clust_fitted = fitter.flatfork_galpy_final_fitting(table, clust_to_fit, iterations, time_to_integrate,
                                                  number_of_steps, try_load=True, graph=False,
                                                  load_fit=False, try_save=False, debug_graph=None)

        # Save it
        with open(windows_directories.orbitsfitdir + "\\" + arrayinfo[0] + "_" +
                  arrayinfo[1] + "_flatfork_fitted_orbit_" + str(clust_to_fit) + ".txt", "wb") as f:
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

def flatfork_do_orbifit_maindata(parameters):

    """
    :param parameters: arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps
    :return:
    """
    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # Grab the clustering
    with open(windows_directories.clusterdir + "\\" + "flatfork_fullgroup_cluster.txt", 'rb') as f:
        clustered = pickle.load(file=f)

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit. Avoid saving until you have finalised the maindata clustering.
        clust_fitted = fitter.flatfork_galpy_fitting_nomemb(table, clustered, clust_to_fit, iterations, time_to_integrate,
                                                   number_of_steps, try_load=False, graph=False,
                                                   load_fit=False, try_save=False,
                                                   extra_text="flatfork_orbifit_maindata_run")

        # Save it
        with open(windows_directories.orbitsfitdir + "\\" + arrayinfo[0] + "_" +
                  arrayinfo[1] + "_flatfork_fitted_orbit_maindata_" + str(clust_to_fit) + ".txt", "wb") as f:
            pickle.dump(obj=clust_fitted, file=f)

def finetune_do_orbifit_maindata(parameters):

    """
    :param parameters: arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps
    :return:
    """
    import energistics

    # Set up parameters and clusters_to_cluster
    arrayinfo, table, clusters_to_cluster, iterations, time_to_integrate, number_of_steps = parameters

    # Set up the fitter
    fitter = energistics.orbifitter()

    # Grab the clustering
    with open(windows_directories.clusterdir_fine + "\\" + "finetune_fullgroup.cluster.txt", 'rb') as f:
        clustered = pickle.load(file=f)

    # For each cluster to cluster
    for clust_to_fit in clusters_to_cluster:

        # Run fit. Avoid saving until you have finalised the maindata clustering.
        clust_fitted = fitter.finetune_galpy_fitting_nomemb(table, clustered, clust_to_fit, iterations, time_to_integrate,
                                                   number_of_steps, try_load=False, graph=False,
                                                   load_fit=False, try_save=False,
                                                   extra_text="finetune_orbifit_maindata_run")

        # Save it
        with open(windows_directories.orbitsfitdir + "\\" + arrayinfo[0] + "_" +
                  arrayinfo[1] + "_finetune_fitted_orbit_maindata_" + str(clust_to_fit) + ".txt", "wb") as f:
            pickle.dump(obj=clust_fitted, file=f)

def maindata_do_orbistatistics(clust_to_fit):

    """

    Run orbistatistics for a given cluster with maindata

    :param clust_to_fit: Which cluster to run orbistatistics for
    :return: stats (see windows_orbifit_maindata)

    """

    # The fitter
    from energistics import orbifitter
    orbifit = orbifitter()

    # Specify savedir/savename and make path, for images
    savedir = os.path.join(os.path.join(windows_directories.imgdir,"orbit_fitting_variables_maindata"),str(clust_to_fit) + "_orbifit_maindata_run")
    save_unique = str(clust_to_fit)
    try:
        os.mkdir(savedir)
    except:
        pass

    # Get orbits for this cluster
    orbits = []

    # In n_carlo
    group = ascii_info.fullgroup
    for saveid in ascii_info.flatfork_orbifit_maindata_saveids:

        # Load it
        with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                  saveid + "_fitted_orbit_maindata_" + str(clust_to_fit) + ".txt", "rb") as f:
            clust_fitted = pickle.load(file=f)
            orbits.append(clust_fitted)

    # Run stats
    stats = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)
    ees, meanee, stdee, \
    periggs, meanpg, stdpg, \
    EEs, meanEE, stdEE, \
    Lzs, meanLz, stdLz, \
    apoggs, meanapog, stdapog = stats

    # Produce various plots
    twod_graph = graphutils.twod_graph()
    twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
    twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
    twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_energy")
    twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_Lz")
    twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Return the stats
    return stats

def flatfork_maindata_do_orbistatistics(clust_to_fit):

    """

    Run orbistatistics for a given cluster with maindata

    :param clust_to_fit: Which cluster to run orbistatistics for
    :return: stats (see windows_orbifit_maindata)

    """

    # The fitter
    from energistics import orbifitter
    orbifit = orbifitter()

    # Specify savedir/savename and make path, for images
    savedir = os.path.join(os.path.join(windows_directories.imgdir,"flatfork_orbit_fitting_variables_maindata"),str(clust_to_fit) + "_orbifit_maindata_run")
    save_unique = str(clust_to_fit)
    try:
        os.mkdir(savedir)
    except:
        pass

    # Get orbits for this cluster
    orbits = []

    # In n_carlo
    group = ascii_info.fullgroup
    for saveid in ascii_info.flatfork_orbifit_maindata_saveids:

        # Load it
        with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                  saveid + "_flatfork_fitted_orbit_maindata_" + str(clust_to_fit) + ".txt", "rb") as f:
            clust_fitted = pickle.load(file=f)
            orbits.append(clust_fitted)

    # Run stats
    stats = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)
    ees, meanee, stdee, \
    periggs, meanpg, stdpg, \
    EEs, meanEE, stdEE, \
    Lzs, meanLz, stdLz, \
    apoggs, meanapog, stdapog = stats

    # Produce various plots
    twod_graph = graphutils.twod_graph()
    twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
    twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
    twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_energy")
    twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_Lz")
    twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Return the stats
    return stats

def finetune_maindata_do_orbistatistics(clust_to_fit):

    """

    Run orbistatistics for a given cluster with maindata

    :param clust_to_fit: Which cluster to run orbistatistics for
    :return: stats (see windows_orbifit_maindata)

    """

    # The fitter
    from energistics import orbifitter
    orbifit = orbifitter()

    # Specify savedir/savename and make path, for images
    savedir = os.path.join(os.path.join(windows_directories.imgdir,"finetune_orbit_fitting_variables_maindata"),str(clust_to_fit) + "_orbifit_maindata_run")
    save_unique = str(clust_to_fit)
    try:
        os.mkdir(savedir)
    except:
        pass

    # Get orbits for this cluster
    orbits = []

    # In n_carlo
    group = ascii_info.fullgroup
    for saveid in ascii_info.finetune_orbifit_maindata_saveids:

        # Load it
        with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                  saveid + "_finetune_fitted_orbit_maindata_" + str(clust_to_fit) + ".txt", "rb") as f:
            clust_fitted = pickle.load(file=f)
            orbits.append(clust_fitted)

    # Run stats
    stats = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)
    ees, meanee, stdee, \
    periggs, meanpg, stdpg, \
    EEs, meanEE, stdEE, \
    Lzs, meanLz, stdLz, \
    apoggs, meanapog, stdapog = stats

    # Produce various plots
    twod_graph = graphutils.twod_graph()
    twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
    twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
    twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_energy")
    twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_Lz")
    twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Return the stats
    return stats


def do_orbistatistics(clust_to_fit):
    """

    Run orbistatistics for a given cluster with memberpercent greatcircle data

    :param clust_to_fit: Which cluster to run orbistatistics for
    :return: stats (see windows_orbifit_arti..._greatcount)

    """

    # The fitter
    from energistics import orbifitter
    orbifit = orbifitter()

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
    group = ascii_info.fullgroup
    saveids = ascii_info.orbifit_saveids

    for saveid in saveids:

        # Load it
        with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                  saveid + "_fitted_orbit_" + str(clust_to_fit) + ".txt", "rb") as f:
            clust_fitted = pickle.load(file=f)
            orbits.append(clust_fitted)

    # Run stats
    stats = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)
    ees, meanee, stdee, \
    periggs, meanpg, stdpg, \
    EEs, meanEE, stdEE, \
    Lzs, meanLz, stdLz, \
    apoggs, meanapog, stdapog = stats

    # Produce various plots
    twod_graph = graphutils.twod_graph()
    twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
    twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
    twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$',
                          savedir + "\\" + str(clust_to_fit) + "_energy")
    twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$',
                          savedir + "\\" + str(clust_to_fit) + "_Lz")
    twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Return the stats
    return stats

def finetune_do_orbistatistics(clust_to_fit):
    """

    Run orbistatistics for a given cluster with memberpercent greatcircle data

    :param clust_to_fit: Which cluster to run orbistatistics for
    :return: stats (see windows_orbifit_arti..._greatcount)

    """

    # The fitter
    from energistics import orbifitter
    orbifit = orbifitter()

    # Specify savedir/savename and make path, for images
    try:
        os.mkdir(windows_directories.imgdir + "\\finetune_orbit_fitting_variables")
    except:
        pass
    savedir = windows_directories.imgdir + "\\finetune_orbit_fitting_variables" + "\\" + str(clust_to_fit)
    save_unique = str(clust_to_fit)
    try:
        os.mkdir(savedir)
    except:
        pass

    # Get orbits for this cluster
    orbits = []

    # In n_carlo
    group = ascii_info.fullgroup
    saveids = ascii_info.orbifit_saveids

    for saveid in saveids:

        # Load it
        with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                  saveid + "_finetune_fitted_orbit_" + str(clust_to_fit) + ".txt", "rb") as f:
            clust_fitted = pickle.load(file=f)
            orbits.append(clust_fitted)

    # Run stats
    stats = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)
    ees, meanee, stdee, \
    periggs, meanpg, stdpg, \
    EEs, meanEE, stdEE, \
    Lzs, meanLz, stdLz, \
    apoggs, meanapog, stdapog = stats

    # Produce various plots
    twod_graph = graphutils.twod_graph()
    twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
    twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
    twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$',
                          savedir + "\\" + str(clust_to_fit) + "_energy")
    twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$',
                          savedir + "\\" + str(clust_to_fit) + "_Lz")
    twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Return the stats
    return stats


def flatfork_do_orbistatistics(clust_to_fit):
    """

    Run orbistatistics for a given cluster with memberpercent greatcircle data

    :param clust_to_fit: Which cluster to run orbistatistics for
    :return: stats (see windows_orbifit_arti..._greatcount)

    """

    # The fitter
    from energistics import orbifitter
    orbifit = orbifitter()

    # Specify savedir/savename and make path, for images
    savedir = windows_directories.imgdir + "\\flatfork_orbit_fitting_variables" + "\\" + str(clust_to_fit)
    save_unique = str(clust_to_fit)
    try:
        os.mkdir(savedir)
    except:
        pass

    # Get orbits for this cluster
    orbits = []

    # In n_carlo
    group = ascii_info.fullgroup
    saveids = ascii_info.orbifit_saveids

    for saveid in saveids:

        # Load it
        with open(windows_directories.orbitsfitdir + "\\" + group + "_" +
                  saveid + "_flatfork_fitted_orbit_" + str(clust_to_fit) + ".txt", "rb") as f:
            clust_fitted = pickle.load(file=f)
            orbits.append(clust_fitted)

    # Run stats
    stats = orbifit.orbistatistics(orbits, 0.3e9, 1000, savedir, save_unique)
    ees, meanee, stdee, \
    periggs, meanpg, stdpg, \
    EEs, meanEE, stdEE, \
    Lzs, meanLz, stdLz, \
    apoggs, meanapog, stdapog = stats

    # Produce various plots
    twod_graph = graphutils.twod_graph()
    twod_graph.hist_fancy(ees, 10, "e", r'$\rho$', savedir + "\\" + str(clust_to_fit) + "_eccentricity")
    twod_graph.hist_fancy(periggs, 10, "perigalacticon / kpc", r'$\rho$', savedir + "\\" + "perigalacticon")
    twod_graph.hist_fancy(EEs, 10, "E / " + r'$\textrm{km}^2\textrm{s}^{-2}$', r'$\rho$',
                          savedir + "\\" + str(clust_to_fit) + "_energy")
    twod_graph.hist_fancy(Lzs, 10, "Lz /" + r'$\textrm{kpc}\textrm{kms}^{-1}$', r'$\rho$',
                          savedir + "\\" + str(clust_to_fit) + "_Lz")
    twod_graph.hist_fancy(apoggs, 10, "apogalacticon / kpc", r'$\rho$', savedir + "\\" + "apogalacticon")

    # Return the stats
    return stats

def flatfork_greatcircle_optimize(parameterss):

    real_dist, which_try, table_included, clustered_included, pole_thetas, pole_phis, resolution = parameterss

    greatfit = galcentricutils.greatfit()

    # Little catching specification for specific clusters (to allow specification on real_dist per-cluster.)
    real_dist_local = real_dist

    # In this case, just do least-squares minimization
    if real_dist_local == True:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmin(least_squares)
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])


    # In this case, maximize the number of stars within real_dist of a greatcircle
    else:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmax(
            least_squares)  # argmax, to mamimize the membership of these cluster members (using a real_dist- see greatfit.least_squares())
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])

    # Take best_circle and obtain the width of the GCC (i.e. get the maximum distance, on the sky, from the GCC.
    std_dist, max_dist = greatfit.deviation_from_gc(table_included['theta'], table_included['phi'], best_circle[0],
                                                    best_circle[1], resolution)

    pole_stdev = np.rad2deg(std_dist)
    pole_maxdist = np.rad2deg(max_dist)

    # Go forth and plot the greatcircle against the data, in ra/dec.
    try:
        graphutils.spec_graph().clust_thetaphi(table=table_included, clustering=clustered_included,
                                               cluster_id=which_try,
                                               vasiliev=False, savedexdir="flatfork_greatcount_beforeMC_" +
                                                                          str(which_try) + "_thetaphi_greatcircle",
                                               gcc_thetaphis=greatfit.gcc_gen(1000, *best_circle),
                                               flatfork=True)
        import matplotlib.pyplot as plt
        plt.close()
        plt.clf()
    except Exception as e:
        print(e)
        pass

    # Return
    return best_circle, pole_stdev, pole_maxdist

def finetune_greatcircle_optimize(parameterss):

    real_dist, which_try, table_included, clustered_included, pole_thetas, pole_phis, resolution = parameterss

    greatfit = galcentricutils.greatfit()

    # Little catching specification for specific clusters (to allow specification on real_dist per-cluster.)
    real_dist_local = real_dist

    # In this case, just do least-squares minimization
    if real_dist_local == True:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmin(least_squares)
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])


    # In this case, maximize the number of stars within real_dist of a greatcircle
    else:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmax(
            least_squares)  # argmax, to mamimize the membership of these cluster members (using a real_dist- see greatfit.least_squares())
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])

    # Take best_circle and obtain the width of the GCC (i.e. get the maximum distance, on the sky, from the GCC.
    std_dist, max_dist = greatfit.deviation_from_gc(table_included['theta'], table_included['phi'], best_circle[0],
                                                    best_circle[1], resolution)

    pole_stdev = np.rad2deg(std_dist)
    pole_maxdist = np.rad2deg(max_dist)

    # Go forth and plot the greatcircle against the data, in ra/dec.
    try:
        graphutils.spec_graph().clust_thetaphi(table=table_included, clustering=clustered_included,
                                               cluster_id=which_try,
                                               vasiliev=False, savedexdir="finetune_greatcount_beforeMC_" +
                                                                          str(which_try) + "_thetaphi_greatcircle",
                                               gcc_thetaphis=greatfit.gcc_gen(1000, *best_circle),
                                               flatfork="finetune")
        import matplotlib.pyplot as plt
        plt.close()
        plt.clf()
    except Exception as e:
        print(e)
        pass

    # Return
    return best_circle, pole_stdev, pole_maxdist

def flatfork_greatcircle_optimize_memberpercent(parameterss):

    real_dist, which_try, table_included, clustered_included, pole_thetas, pole_phis, resolution = parameterss

    greatfit = galcentricutils.greatfit()

    # Little catching specification for specific clusters (to allow specification on real_dist per-cluster.)
    real_dist_local = real_dist

    # In this case, just do least-squares minimization
    if real_dist_local == True:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmin(least_squares)
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])


    # In this case, maximize the number of stars within real_dist of a greatcircle
    else:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmax(
            least_squares)  # argmax, to mamimize the membership of these cluster members (using a real_dist- see greatfit.least_squares())
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])

    # Take best_circle and obtain the width of the GCC (i.e. get the maximum distance, on the sky, from the GCC.
    std_dist, max_dist = greatfit.deviation_from_gc(table_included['theta'], table_included['phi'], best_circle[0],
                                                    best_circle[1], resolution)

    pole_stdev = np.rad2deg(std_dist)
    pole_maxdist = np.rad2deg(max_dist)

    # Go forth and plot the greatcircle against the data, in ra/dec.
    try:
        graphutils.spec_graph().clust_thetaphi(table=table_included, clustering=clustered_included,
                                               cluster_id=which_try,
                                               vasiliev=False, savedexdir="flatfork_greatcount_" +
                                                                          str(which_try) + "_thetaphi_greatcircle",
                                               gcc_thetaphis=greatfit.gcc_gen(1000, *best_circle),
                                               flatfork=True)
        import matplotlib.pyplot as plt
        plt.close()
        plt.clf()
    except Exception as e:
        print(e)
        pass

    # Return
    return best_circle, pole_stdev, pole_maxdist

def finetune_greatcircle_optimize_memberpercent(parameterss):

    real_dist, which_try, table_included, clustered_included, pole_thetas, pole_phis, resolution = parameterss

    greatfit = galcentricutils.greatfit()

    # Little catching specification for specific clusters (to allow specification on real_dist per-cluster.)
    real_dist_local = real_dist

    # In this case, just do least-squares minimization
    if real_dist_local == True:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmin(least_squares)
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])


    # In this case, maximize the number of stars within real_dist of a greatcircle
    else:

        least_squares = []

        for theta_pole, phi_pole in zip(pole_thetas, pole_phis):
            leastsq = greatfit.least_squares(table_included['theta'], table_included['phi'], theta_pole, phi_pole,
                                             resolution, real_dist_local)
            least_squares.append(leastsq)

        best_circle = np.argmax(
            least_squares)  # argmax, to mamimize the membership of these cluster members (using a real_dist- see greatfit.least_squares())
        best_circle = np.array([pole_thetas[best_circle], pole_phis[best_circle]])

    # Take best_circle and obtain the width of the GCC (i.e. get the maximum distance, on the sky, from the GCC.
    std_dist, max_dist = greatfit.deviation_from_gc(table_included['theta'], table_included['phi'], best_circle[0],
                                                    best_circle[1], resolution)

    pole_stdev = np.rad2deg(std_dist)
    pole_maxdist = np.rad2deg(max_dist)

    # Go forth and plot the greatcircle against the data, in ra/dec.
    try:
        graphutils.spec_graph().clust_thetaphi(table=table_included, clustering=clustered_included,
                                               cluster_id=which_try,
                                               vasiliev=False, savedexdir="flatfork_greatcount_" +
                                                                          str(which_try) + "_thetaphi_greatcircle",
                                               gcc_thetaphis=greatfit.gcc_gen(1000, *best_circle),
                                               flatfork="finetune")
        import matplotlib.pyplot as plt
        plt.close()
        plt.clf()
    except Exception as e:
        print(e)
        pass

    # Return
    return best_circle, pole_stdev, pole_maxdist


# Factor by which to multiply max_clust_size of the GSE (to prevent central clumps being caught by GSE.)
massive_factor = 1.1
def flatfork_do_finetune(args):

    from galcentricutils import compclust, cluster3d
    cc = compclust()
    c3d = cluster3d()

    # Args
    cut_table, clust_to_finetune, finetune_args = args

    # Args
    minsamples_range, \
    min_clust_size, \
    trials_per_samples, \
    minimum_trial_score, \
    minimum_trial_analytically = finetune_args

    # Isolate the gcc parameters
    cut_clustered = cut_table['prelim_clust']

    # Determine size of GSE (or the largest cluster in the original sample- GSE default.)
    is_gse = np.array([True if d == ascii_info.flatfork_GSE_ID else False for d in cut_clustered])
    is_gse = np.sum(is_gse) * (massive_factor)
    cut_array = np.array([cut_table['Lx'], cut_table['Ly'], cut_table['Lz']]).T
    max_clust_size = is_gse
    best_samples, passfail, \
    numclusts, scores, sizedifs = c3d.minsamples_finetune(cut_array, cut_clustered,
                                                          clust_to_finetune, minsamples_range, min_clust_size,
                                                          trials_per_samples, minimum_trial_score,
                                                          minimum_trial_analytically, max_clust_size)

    # Cluster it accordingly (with the new samples- for debugging/etc)
    from hdbscan import HDBSCAN
    trial_clustered = HDBSCAN(min_cluster_size=int(min_clust_size),
                              max_cluster_size=is_gse,
                              min_samples=int(best_samples),
                              algorithm='best',
                              metric='l2').fit_predict(cut_array)

    # Graphing (for debug.)
    import graphutils
    finetune_debugdir = "flatfork_finetune"
    finetune_path_gcc = os.path.join(finetune_debugdir, str(clust_to_finetune) + "_gcc")
    finetune_path_newsamples = os.path.join(finetune_debugdir, str(clust_to_finetune) + "_gcc_newsamples")
    graphutils.threed_graph().kmeans_L_array(cut_array, cut_clustered,
                                             finetune_path_gcc, browser=False, outliers=True)
    trial_clustered = cc.compclust_multilabel_julia(cut_clustered, trial_clustered, 30)
    graphutils.threed_graph().kmeans_L_array(cut_array, trial_clustered,
                                             finetune_path_newsamples, browser=False, outliers=True)

    # Evaluate the sums for the initial GCC and the final GCC with newsamples
    inisum = len(np.where(cut_clustered==clust_to_finetune)[0])
    finsum = len(np.where(trial_clustered==clust_to_finetune)[0])

    # Return
    print(best_samples, inisum, finsum, clust_to_finetune)
    return best_samples, inisum, finsum

def finetune_do_finetune(args):

    from galcentricutils import compclust, cluster3d
    cc = compclust()
    c3d = cluster3d()

    # Args
    cut_table, clust_to_finetune, finetune_args = args

    # Args
    minsamples_range, \
    min_clust_size, \
    trials_per_samples, \
    minimum_trial_score, \
    minimum_trial_analytically = finetune_args

    # Isolate the gcc parameters
    cut_clustered = cut_table['prelim_clust']

    # Determine size of GSE (or the largest cluster in the original sample- GSE default.)
    is_gse = np.array([True if d == ascii_info.finetune_GSE_ID else False for d in cut_clustered])
    is_gse = np.sum(is_gse) * (massive_factor)
    cut_array = np.array([cut_table['Lx'], cut_table['Ly'], cut_table['Lz']]).T
    max_clust_size = is_gse
    best_samples, passfail, \
    numclusts, scores, sizedifs = c3d.minsamples_finetune(cut_array, cut_clustered,
                                                          clust_to_finetune, minsamples_range, min_clust_size,
                                                          trials_per_samples, minimum_trial_score,
                                                          minimum_trial_analytically, max_clust_size)

    # Cluster it accordingly (with the new samples- for debugging/etc)
    from hdbscan import HDBSCAN
    trial_clustered = HDBSCAN(min_cluster_size=int(min_clust_size),
                              max_cluster_size=is_gse,
                              min_samples=int(best_samples),
                              algorithm='best',
                              metric='l2').fit_predict(cut_array)

    # Graphing (for debug.)
    import graphutils
    finetune_debugdir = "finetune_finetune"
    finetune_path_gcc = os.path.join(finetune_debugdir, str(clust_to_finetune) + "_gcc")
    finetune_path_newsamples = os.path.join(finetune_debugdir, str(clust_to_finetune) + "_gcc_newsamples")
    graphutils.threed_graph().kmeans_L_array(cut_array, cut_clustered,
                                             finetune_path_gcc, browser=False, outliers=True)
    trial_clustered = cc.compclust_multilabel_julia(cut_clustered, trial_clustered, 30)
    graphutils.threed_graph().kmeans_L_array(cut_array, trial_clustered,
                                             finetune_path_newsamples, browser=False, outliers=True)
    import matplotlib.pyplot as plt
    plt.close()
    plt.clf()

    # Orbifit it with the new clustered too! Only for smaller ones (less than a thousand stars.)
    where_to_finetune = np.where(trial_clustered==clust_to_finetune)[0]
    if len(where_to_finetune) < 1000 and clust_to_finetune in trial_clustered:

        cut_table['prelim_clust'] = trial_clustered

        orbifit = energistics.orbifitter().finetune_galpy_fitting_nomemb(cut_table,
                                                                         trial_clustered,
                                                                         clust_to_finetune,
                                                                         2000,
                                                                         0.6e9,
                                                                         1000,
                                                                         False,
                                                                         False,
                                                                         False,
                                                                         False,
                                                                         extra_text="for_lbplots_finetune_debug")



        # Forward Integral
        import copy
        from astropy import units as u
        forward = copy.deepcopy(orbifit)
        backward = copy.deepcopy(orbifit)
        forward.integrate((np.linspace(0, 1e9, 2000) * u.yr), energistics.orbigistics().pot)
        backward.integrate((np.linspace(0, -1e9, 2000) * u.yr), energistics.orbigistics().pot)
        llsf, bbsf, llsb, bbsb = forward.ll((np.linspace(0, 0.7e9, 2000) * u.yr)).value, \
                                 forward.bb((np.linspace(0, 0.7e9, 2000) * u.yr)).value, \
                                 backward.ll((np.linspace(0, -0.15e9, 2000) * u.yr)).value, \
                                 backward.bb((np.linspace(0, -0.15e9, 2000) * u.yr)).value
        lls, bbs = np.concatenate([llsf, llsb]), np.concatenate([bbsf, bbsb])
        lls = [d - 360 if d > 180 else d for d in lls]

        specgrapher = graphutils.spec_graph()
        fig, axs = specgrapher.lb_orbits(cut_table[[True if d == clust_to_finetune
                                                else False
                                                for d in trial_clustered]],
                                         0.3e9, [-180, 180],
                                         [-90, 90], None,
                                         line=False, points=4000)
        axs.scatter(lls, bbs, color='lime', marker='o', s=3)
        # Misc axis things
        axs.set_facecolor("k")
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color='white')
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color='white')
        axs.grid(color="white")

        # Savedir
        savepath = os.path.join(os.path.join(windows_directories.imgdir, finetune_debugdir), str(clust_to_finetune) + ".png")
        import matplotlib.pyplot as plt
        plt.savefig(savepath, dpi=300, transparent=False)
        plt.close()

    # Evaluate the sums for the initial GCC and the final GCC with newsamples
    inisum = len(np.where(cut_clustered==clust_to_finetune)[0])
    finsum = len(np.where(trial_clustered==clust_to_finetune)[0])

    # Return
    print("finetune_finetune", best_samples, inisum, finsum, clust_to_finetune)
    return best_samples, inisum, finsum


# Minimum entropy fit
def do_minentropy(parameters):

    table, clustering, streams_todo, ranges, nmonte = parameters


    # Set up an empty table to hold detail for each cluster (this assumes we do this per cluster.)
    pointtable = Table(names=["clust", "M_d", "M_nfw", "c_nfw"], data=np.array([streams_todo,
                                                                                np.zeros_like(streams_todo),
                                                                                np.zeros_like(streams_todo),
                                                                                np.zeros_like(streams_todo)]).T)
    pointtable['clust'] = streams_todo

    # Generate entropies and save.
    points_saved = []
    for clusts_to_fit in streams_todo:
        clusts_to_fit = [clusts_to_fit]
        # Grab and produce required data lists
        x,y,z,vx,vy,vz = [[] for d in range(6)]
        for clust_to_fit in clusts_to_fit:
            data_to_fit = table[[True if d == clust_to_fit else False for d in clustering]]
            xx, yy, zz = data_to_fit['x'], data_to_fit['y'], data_to_fit['z']
            vxx, vyy, vzz = data_to_fit['vx'], data_to_fit['vy'], data_to_fit['vz']
            x.append(xx)
            y.append(yy)
            z.append(zz)
            vx.append(vxx)
            vy.append(vyy)
            vz.append(vzz)

        # Set up entrofitter
        entrofit = energistics.entro_fitter(x,y,z,
                                            vx,vy,vz,
                                            [ranges[0][0]*M_d, ranges[0][1]*M_d], [ranges[1][0]*M_nfw, ranges[1][1]*M_nfw], [ranges[2][0]*c_nfw, ranges[2][1]*c_nfw],
                                            [M_b, a_b, M_d, a_d, b_d, M_nfw, a_nfw, c_nfw],
                                            nmonte,
                                            0.2, 100)

        # Set up points
        points = entrofit.genpoints()
        entropies = entrofit.list_entropy(points)

        # Get the minimum point
        best_point = points[np.argmin(entropies)]/np.array([M_d, M_nfw, c_nfw])
        points_saved.append(best_point)

    points_saved = np.array(points_saved).T
    pointtable['M_d'], pointtable['M_nfw'], pointtable['c_nfw'] = points_saved
    return pointtable

# The above, but combines all the provided clusters into one (i.e. multiple stream combination).
# Here streamstodo is a list of lists, i.e. [[1,2,3,4], [4,5,6,7]], with entropy min over sum of each list
def do_seventropy(parameters):

    table, clustering, streams_todo, ranges, nmonte = parameters

    # Generate entropies and save.
    points_saved = []

    for clusts_to_fit in streams_todo:

        # Grab and produce required data lists
        x,y,z,vx,vy,vz = [],[],[],[],[],[]
        for clust_to_fit in clusts_to_fit:
            data_to_fit = table[[True if d == clust_to_fit else False for d in clustering]]
            xx, yy, zz = data_to_fit['x'], data_to_fit['y'], data_to_fit['z']
            vxx, vyy, vzz = data_to_fit['vx'], data_to_fit['vy'], data_to_fit['vz']
            x.append(xx)
            y.append(yy)
            z.append(zz)
            vx.append(vxx)
            vy.append(vyy)
            vz.append(vzz)

        # Set up entrofitter
        entrofit = energistics.entro_fitter(x,y,z,
                                            vx,vy,vz,
                                            [ranges[0][0]*M_d, ranges[0][1]*M_d], [ranges[1][0]*M_nfw, ranges[1][1]*M_nfw], [ranges[2][0]*c_nfw, ranges[2][1]*c_nfw],
                                            [M_b, a_b, M_d, a_d, b_d, M_nfw, a_nfw, c_nfw],
                                            nmonte,
                                            0.2, 100)

        # Set up points
        points = entrofit.genpoints()
        entropies = entrofit.list_entropy(points)

        # Get the minimum point
        best_point = points[np.argmin(entropies)]/np.array([M_d, M_nfw, c_nfw])
        points_saved.append(best_point)

    best_point = np.array(points_saved)[0]

    return best_point

