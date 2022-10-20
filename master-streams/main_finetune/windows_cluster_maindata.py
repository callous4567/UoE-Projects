import os
import pickle
import numpy as np
import ascii_info
import galcentricutils
import hdfutils
import windows_directories
from hdbscan import flat
import graphutils
from matplotlib import pyplot as plt
from energistics import orbifitter
from windows_directories import clusterer_flat_path, clusterer_flat_labels

plotfit = True
iterations, time_to_integrate, number_of_steps, try_load, graph = 250, 0.3e9, 250, False, True
load_fit, try_save = False, False
max_size = 250 # of cluster
vasiliev = True
runfit = True
fine_nclust = ascii_info.fine_nclust # from the MSc project we ascertained this is a good number based on results
writer = hdfutils.hdf5_writer(windows_directories.datadir,
                              ascii_info.finetune_asciiname)

if runfit == True:

    # Clustering parameters/etc
    arrayinfominpars = []
    group = ascii_info.fullgroup
    minpar = ascii_info.fulldata_minpar

    # Set up data and run the clustering (
    panda = writer.read_df(ascii_info.fullgroup,ascii_info.panda_raw)

    """
    # Testing with some slightly noised data (we are after all, only trying to estimate the locus of the cluster.)
    rng = numpy.random.default_rng()
    for vec_L in panda['vec_L']:
        vec_L += rng.normal(loc=0, scale=100, size=3) # REDUNDANT! But does improve results slightly. 
    """
    data = np.array(panda['vec_L'])
    data = list(data)

    # Set up the flat hdbscan run. 25 works well (See Finetune- 25/26.) This is more predictable to produce a mean.
    from hdbscan import flat
    clusterer = flat.HDBSCAN_flat(X=data,
                                  n_clusters=fine_nclust,
                                  min_cluster_size=minpar[0],
                                  min_samples=minpar[1],
                                  metric='l2',
                                  algorithm='best',
                                  prediction_data=False)

    clustered = clusterer.labels_
    with open(windows_directories.clusterdir_fine + "\\" + "finetune_fullgroup.cluster.txt", 'wb') as f:
        pickle.dump(obj=clustered, file=f)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found.
    numclust = galcentricutils.compclust().nclust_get(clustered)
    savedexdir = "clustered_finetune\\" + "finetune_fullgroup_clustered_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=True,  outliers=False)

    # Generate a HTML for this clustering under imgdir\\clustered and embed the number of clusters found with the outliers
    numclust = galcentricutils.compclust().nclust_get(clustered)
    savedexdir = "clustered_finetune\\" + "finetune_fullgroup_clustered_outliers_" + str(numclust) + "_clusters_found"
    data = np.array(data)
    graphutils.threed_graph().kmeans_L_array(data, clustered, savedexdir, browser=True,  outliers=True)

    # Generate a table of membership fractions/etc
    fracs = [0 for d in np.arange(-1,fine_nclust,1)]
    for num,cluster in enumerate(np.arange(-1,fine_nclust,1)):
        for star in clustered:
            if star == cluster:
                fracs[num] += 1
    fracs = np.array(fracs, dtype=float)
    fracs /= len(clustered)
    fracs_array = np.array([np.arange(-1, fine_nclust, 1), fracs]).T
    writer.write(ascii_info.fullgroup,"fine_fractions", fracs_array)

    # Load Data
    table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)

    # Save preliminary clustering to the table
    table['prelim_clust'] = clustered
    writer.write_table(ascii_info.fullgroup,ascii_info.fullset,table)

table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)

clustered = np.array(table['prelim_clust'], int)
numclust = galcentricutils.compclust().nclust_get(clustered)

if vasiliev == True:

    # Test out the vasiliev graphing for this set...
    table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)

    for clust_id in range(numclust):
        #if clust_id == 13:
        #    graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
        #                                        savedexdir="finetune_" + str(clust_id) + "_clusttest_lb", lb=True, vasiliev=True)
        #    graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
        #                                        savedexdir="finetune_" + str(clust_id) + "_clusttest_ra", lb=False, vasiliev=True)
        #else:
        try:
            graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                savedexdir="finetune_" + str(clust_id) + "_clusttest_lb",
                                                lb=True,
                                                vasiliev=False,
                                                flatfork="finetune")
            graphutils.spec_graph().clust_radec(table, clustered, cluster_id=clust_id,
                                                savedexdir="finetune_" + str(clust_id) + "_clusttest_ra",
                                                lb=False,
                                                vasiliev=False,
                                                flatfork="finetune")
        except:
            pass

        # Also generate regular lb plots
        #savepath = windows_directories.imgdir + "\\" + "vasiliev"+ "\\" + str(clust_id) + "_clusttest_lbplot" + ".png"
        #graphutils.twod_graph().lbplot(table[[True if d == clust_id else False for d in clustered]], savepath,
        #                               negpi=True)

if plotfit == True:

    iterations, time_to_integrate, number_of_steps, try_load, graph = 250, 0.3e9, 250, False, False
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
                                               number_of_steps, False, True, False, False)
            try:
                os.mkdir(windows_directories.imgdir + "\\finetune_orbit_fitting_variables_maindata")
            except:
                pass

            savedir = windows_directories.imgdir + "\\finetune_orbit_fitting_variables_maindata" + "\\" + str(cluster)

            try:
                os.mkdir(savedir)
            except:
                pass

            specgrapher.orbiplot(table, clustered, cluster, fit, 0.3e9, 250, savedir, "practice_plot")
