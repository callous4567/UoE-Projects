import os
import pickle
import numpy as np
import ascii_info_new
import galcentricutils_new
import hdfutils
import windows_directories_new
from hdbscan import flat
import graphutils_new
from matplotlib import pyplot as plt
from windows_directories_new import clusterer_flat_path, clusterer_flat_labels

plotfit = True
iterations, time_to_integrate, number_of_steps, try_load, graph = 250, 0.3e9, 250, False, True
load_fit, try_save = False, False
max_size = 250 # of cluster
vasiliev = True
runflat = False
flat_nclust = 24 # from the MSc project we ascertained this is a good number based on results

if runflat == True:

    # Clustering parameters/etc
    arrayinfominpars = []
    group = ascii_info_new.fullgroup
    minpar = ascii_info_new.fulldata_minpar

    # Set up data and run the clustering
    panda = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_df(ascii_info_new.fullgroup,
                                                                        ascii_info_new.panda_raw)

    """
    # Set up slightly-noised data (noise to ensure any irregularities in data distribution are removed)
    rng = numpy.random.default_rng()
    for vec_L in panda['vec_L']:
        vec_L += rng.normal(loc=0, scale=0, size=3)
    """
    data = np.array(panda['vec_L'])
    data = list(data)


    # Set up the flat hdbscan run
    clusterer = flat.HDBSCAN_flat(X=data,
                                  n_clusters=flat_nclust,
                                  min_cluster_size=minpar[0],
                                  min_samples=minpar[1],
                                  metric='l2',
                                  algorithm='best',
                                  prediction_data=True)

    # Save the clusterer
    with open(clusterer_flat_path, 'wb') as f:
        pickle.dump(obj=clusterer, file=f)

    # Save the label for compatibility, too
    with open(clusterer_flat_labels, 'wb') as f:
        pickle.dump(obj=clusterer.labels_, file=f)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
    numclust = galcentricutils_new.compclust().nclust_get(clusterer.labels_)
    savedexdir = "\\clustered\\flatfork_fullgroup_clustered_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils_new.threed_graph().kmeans_L_array(data, clusterer.labels_, savedexdir, browser=True,  outliers=False)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found with the outliers
    numclust = galcentricutils_new.compclust().nclust_get(clusterer.labels_)
    savedexdir = "\\clustered\\flatfork_fullgroup_clustered_outliers_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils_new.threed_graph().kmeans_L_array(data, clusterer.labels_, savedexdir, browser=True,  outliers=True)

    # Generate a table of membership fractions/etc
    fracs = [0 for d in np.arange(-1,np.max(clusterer.labels_) + 1,1)]
    for num,cluster in enumerate(np.arange(-1,np.max(clusterer.labels_) + 1,1)):
        for star in clusterer.labels_:
            if star == cluster:
                fracs[num] += 1
    fracs = np.array(fracs, dtype=float)
    fracs /= len(clusterer.labels_)
    fracs_array = np.array([np.arange(-1, np.max(clusterer.labels_) + 1, 1), fracs]).T
    hdfutils.hdf5_writer(windows_directories_new.datadir,
                         ascii_info_new.flatfork_asciiname).write(ascii_info_new.fullgroup,"flat_fractions", fracs_array)

    # Load Data
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                  ascii_info_new.fullset)

    # Save preliminary clustering to the table
    table['prelim_clust'] = clusterer.labels_
    hdfutils.hdf5_writer(windows_directories_new.datadir,
                         ascii_info_new.flatfork_asciiname).write_table(ascii_info_new.fullgroup,
                                                           ascii_info_new.fullset,
                                                           table)

if plotfit == True:

    # Load Data
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                  ascii_info_new.fullset)
    clustered = np.array(table['prelim_clust'], int)
    numclust = galcentricutils_new.compclust().nclust_get(clustered)

    # Save the clusterer
    with open(clusterer_flat_path, 'rb') as f:
        clusterer = pickle.load(file=f)

    # Set up data and run the clustering
    panda = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_df(ascii_info_new.fullgroup,
                                                                        ascii_info_new.panda_raw)

    """
    # Set up slightly-noised data (noise to ensure any irregularities in data distribution are removed)
    rng = numpy.random.default_rng()
    for vec_L in panda['vec_L']:
        vec_L += rng.normal(loc=0, scale=0, size=3)
    """
    data = np.array(panda['vec_L'])
    data = list(data)


    # Do the HTML plots, too.
    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
    savedexdir = "\\clustered\\flatfork_fullgroup_clustered_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils_new.threed_graph().kmeans_L_array(data, clusterer.labels_, savedexdir, browser=True,  outliers=False)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found with the outliers
    numclust = galcentricutils_new.compclust().nclust_get(clusterer.labels_)
    savedexdir = "\\clustered\\flatfork_fullgroup_clustered_outliers_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils_new.threed_graph().kmeans_L_array(data, clusterer.labels_, savedexdir, browser=True,  outliers=True)

    from energistics_new import orbifitter
    orbifit = orbifitter()
    specgrapher = graphutils_new.spec_graph()
    #clusters_to_try = np.arange(0, numclust, 1)
    maindata_clusters_to_try = []
    for cluster in clustered:
        if cluster not in maindata_clusters_to_try:
            maindata_clusters_to_try.append(cluster)

    for cluster in maindata_clusters_to_try:

        truefalses = [True if d == cluster else False for d in clustered]
        size = np.sum(np.ones(len(truefalses))*np.array(truefalses))

        if size <= max_size:

            print("Cluster size and name of ", size, cluster)

            try:
                plt.close()
            except:
                pass

            fit = orbifit.galpy_fitting_nomemb(table, clustered, cluster, iterations, time_to_integrate,
                                               number_of_steps, try_load, graph, load_fit, try_save, extra_text="flatfork")
            try:
                os.mkdir(windows_directories_new.imgdir + "\\flatfork_orbit_fitting_variables_maindata")
            except:
                pass

            savedir = windows_directories_new.imgdir + "\\" + "flatfork_orbit_fitting_variables_maindata" + "\\" + str(cluster)

            try:
                os.mkdir(savedir)
            except:
                pass

            specgrapher.orbiplot(table, clustered, cluster, fit, 0.3e9, 250, savedir, str(cluster) + "_practice_plot")

            import matplotlib.pyplot as plt
            plt.close()

if vasiliev == True:

    # Load Data
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                  ascii_info_new.fullset)
    clustered = np.array(table['prelim_clust'], int)
    numclust = galcentricutils_new.compclust().nclust_get(clustered)

    for clust_id in range(numclust):

        try:
            graphutils_new.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                savedexdir="flatfork_" + str(clust_id) + "_clusttest_lb",
                                                lb=True, vasiliev=False,
                                                flatfork=True)
            graphutils_new.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                savedexdir="flatfork_" + str(clust_id) + "_clusttest_ra",
                                                lb=False, vasiliev=False,
                                                flatfork=True)
        except:
            pass

