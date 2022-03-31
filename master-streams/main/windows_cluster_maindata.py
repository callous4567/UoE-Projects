import os
import pickle
import numpy as np
from matplotlib import pyplot as plt

import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
import time
from energistics import orbifitter

runfit = True
plotfit = False
vasiliev = False

if runfit == True:
    # Clustering parameters/etc
    arrayinfominpars = []
    group = ascii_info.fullgroup
    minpar = ascii_info.fulldata_minpar

    # Set up data and run the clustering
    panda = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                               ascii_info.panda_raw)
    data = np.array(panda['vec_L'])
    data = list(data)
    clustered = galcentricutils.cluster3d().listhdbs(data, minpar)
    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'wb') as f:
        pickle.dump(obj=clustered, file=f)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
    numclust = galcentricutils.compclust().nclust_get(clustered)
    savedexdir = "\\clustered\\fullgroup_clustered_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=False,  outliers=False)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found with the outliers
    numclust = galcentricutils.compclust().nclust_get(clustered)
    savedexdir = "\\clustered\\fullgroup_clustered_outliers_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=False,  outliers=True)

    # Also do a spatial one
    x,y,z = panda['x'],panda['y'],panda['z']
    pos = np.array([x,y,z]).T
    savedexdir = "\\clustered\\fullgroup_clustered_" + str(numclust) + "_clusters_found_spatial"
    graphutils.threed_graph().xmeans_L_array(pos, clustered, savedexdir, False, False)

    # Generate a table of membership fractions/etc
    fracs = [0 for d in range(np.max(clustered) + 1)]
    for cluster in range(np.max(clustered) + 1):
        for star in clustered:
            if star == cluster:
                fracs[cluster] += 1
    fracs = np.array(fracs)
    fracs /= len(clustered)
    data_array = np.array([np.arange(0, np.max(clustered) + 1, 1), fracs]).T


# Load Data
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustered = pickle.load(file=f)
    numclust = galcentricutils.compclust().nclust_get(clustered)

if vasiliev == True:
    # Test out the vasiliev graphing for this set...
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)
    for clust_id in range(numclust):
        if clust_id == 13:
            graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                savedexdir=str(clust_id) + "_clusttest_lb", lb=True, vasiliev=True)
            graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                savedexdir=str(clust_id) + "_clusttest_ra", lb=False, vasiliev=True)
        else:
            try:
                graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                    savedexdir=str(clust_id) + "_clusttest_lb", lb=True, vasiliev=False)
                graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                    savedexdir=str(clust_id) + "_clusttest_ra", lb=False, vasiliev=False)
            except:
                pass

        # Also generate regular lb plots
        #savepath = windows_directories.imgdir + "\\" + "vasiliev"+ "\\" + str(clust_id) + "_clusttest_lbplot" + ".png"
        #graphutils.twod_graph().lbplot(table[[True if d == clust_id else False for d in clustered]], savepath,
        #                               negpi=True)
if plotfit == True:
    iterations, time_to_integrate, number_of_steps, try_load, graph = 2000, 0.3e9, 1000, True, True
    orbifit = orbifitter()
    specgrapher = graphutils.spec_graph()
    clusters_to_try = np.arange(0, numclust, 1)
    maindata_clusters_to_try = []
    for cluster in clustered:
        if cluster not in maindata_clusters_to_try:
            maindata_clusters_to_try.append(cluster)

    for cluster in maindata_clusters_to_try:
        if cluster in [-1,13,16]:
            pass
        else:
            try:
                plt.close()
            except:
                pass
            fit = orbifit.galpy_fitting_nomemb(table, clustered, cluster, iterations, time_to_integrate,
                                               number_of_steps, True, True, False, False)
            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_maindata" + "\\" + str(cluster)
            specgrapher.orbiplot(table, clustered, cluster, fit, 0.3e9, 1000, savedir, "practice_plot")
