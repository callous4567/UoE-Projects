import os
import pickle
import numpy as np
import numpy.random
from matplotlib import pyplot as plt
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
from energistics import orbifitter

runfit = False
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

    """
    # Testing with some slightly noised data (we are after all, only trying to estimate the locus of the cluster.)
    rng = numpy.random.default_rng()
    for vec_L in panda['vec_L']:
        vec_L += rng.normal(loc=0, scale=100, size=3)
    """

    data = np.array(panda['vec_L'])
    data = list(data)

    # Set up the flat hdbscan run. 25 works well (See Flatfork- 25/26.) This is more predictable to produce a mean.
    from hdbscan import flat
    clusterer = flat.HDBSCAN_flat(X=data,
                                  n_clusters=25,
                                  min_cluster_size=minpar[0],
                                  min_samples=minpar[1],
                                  metric='l2',
                                  algorithm='best',
                                  prediction_data=True)

    clustered = clusterer.labels_
    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'wb') as f:
        pickle.dump(obj=clustered, file=f)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
    numclust = galcentricutils.compclust().nclust_get(clustered)
    savedexdir = "\\clustered\\fullgroup_clustered_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=True,  outliers=False)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found with the outliers
    numclust = galcentricutils.compclust().nclust_get(clustered)
    savedexdir = "\\clustered\\fullgroup_clustered_outliers_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=True,  outliers=True)

    # Also do a spatial one
    x,y,z = panda['x'],panda['y'],panda['z']
    pos = np.array([x,y,z]).T
    savedexdir = "\\clustered\\fullgroup_clustered_" + str(numclust) + "_clusters_found_spatial"
    graphutils.threed_graph().xmeans_L_array(pos, clustered, savedexdir, True, False)

    # Generate a table of membership fractions/etc
    fracs = [0 for d in range(np.max(clustered) + 1)]
    for cluster in range(np.max(clustered) + 1):
        for star in clustered:
            if star == cluster:
                fracs[cluster] += 1
    fracs = np.array(fracs, dtype=float)
    fracs /= len(clustered)
    data_array = np.array([np.arange(0, np.max(clustered) + 1, 1), fracs]).T

    # Load Data
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)

    # Save preliminary clustering to the table
    table['prelim_clust'] = clustered
    hdfutils.hdf5_writer(windows_directories.datadir,
                         ascii_info.asciiname).write_table(ascii_info.fullgroup,
                                                           ascii_info.fullset,
                                                           table)

table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
clustered = np.array(table['prelim_clust'], int)
numclust = galcentricutils.compclust().nclust_get(clustered)

if vasiliev == True:

    # Test out the vasiliev graphing for this set...
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)
    for clust_id in range(numclust):
        #if clust_id == 13:
        #    graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
        #                                        savedexdir=str(clust_id) + "_clusttest_lb", lb=True, vasiliev=True)
        #    graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
        #                                        savedexdir=str(clust_id) + "_clusttest_ra", lb=False, vasiliev=True)
        #else:
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

    iterations, time_to_integrate, number_of_steps, try_load, graph = 250, 0.3e9, 250, True, True
    orbifit = orbifitter()
    specgrapher = graphutils.spec_graph()
    clusters_to_try = np.arange(0, numclust, 1)
    maindata_clusters_to_try = []
    for cluster in clustered:
        if cluster not in maindata_clusters_to_try:
            maindata_clusters_to_try.append(cluster)

    # max size spec
    max_size = 150
    for cluster in maindata_clusters_to_try:

        truefalses = [True if d == cluster else False for d in clustered]
        size = np.sum(np.ones(len(truefalses))*np.array(truefalses))
        if size >= max_size:

            print("Cluster size and name of ", size, cluster)

            try:
                plt.close()
            except:
                pass

            fit = orbifit.galpy_fitting_nomemb(table, clustered, cluster, iterations, time_to_integrate,
                                               number_of_steps, True, True, False, False)
            try:
                os.mkdir(windows_directories.imgdir + "\\orbit_fitting_variables_maindata")
            except:
                pass

            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_maindata" + "\\" + str(cluster)

            try:
                os.mkdir(savedir)
            except:
                pass

            specgrapher.orbiplot(table, clustered, cluster, fit, 0.3e9, 250, savedir, "practice_plot")
